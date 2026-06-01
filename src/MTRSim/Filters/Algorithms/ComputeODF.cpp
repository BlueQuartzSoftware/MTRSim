#include "ComputeODF.hpp"

#include "simplnx/Common/Range.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/Utilities/MaskCompareUtilities.hpp"
#include "simplnx/Utilities/MessageHelper.hpp"
#include "simplnx/Utilities/ParallelDataAlgorithm.hpp"

#include "LibMTRSim/ODFBuilder.hpp"

#include <EbsdLib/Core/EbsdLibConstants.h>
#include <EbsdLib/LaueOps/LaueOps.h>
#include <EbsdLib/Math/Matrix3X3.hpp>
#include <EbsdLib/Orientation/Euler.hpp>
#include <EbsdLib/Orientation/OrientationFwd.hpp>
#include <EbsdLib/Orientation/OrientationMatrix.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <mutex>
#include <numbers>
#include <vector>

using namespace nx::core;

namespace {

// Worker functor invoked by ParallelDataAlgorithm. Each invocation owns a
// private accumulator vector and a private contributing-voxel count. On
// destruction the worker merges those private results into the shared master
// accumulator/count under a mutex.
//
// This pattern avoids any concurrent writes to the master state and matches the
// thread-safety guidance: per-thread (here, per-call) std::vector<double> is
// the only mutable thing touched in the parallel section.
//
// All orientation math goes through EbsdLib's `Euler<double>`,
// `OrientationMatrix<double>`, and `LaueOps` directly. Float32 EBSD inputs
// are promoted to `double` at the function boundary; all subsequent math is
// double precision (matches EbsdLib's internal convention and the project's
// double-precision orientation rule).
class AccumulateWorker {
public:
  AccumulateWorker(const Float32Array &eulers, const Int32Array &phases,
                   const UInt32Array &crystalStructures,
                   const MaskCompareUtilities::MaskCompare *mask,
                   const mtrsim::ODFBuildParams &params,
                   std::vector<double> &masterValues,
                   std::size_t &masterContributingCount,
                   std::size_t &masterTotalDeposits,
                   std::size_t &masterFailureCount, std::mutex &mergeMutex,
                   const std::atomic_bool &shouldCancel,
                   ProgressMessageHelper &progressHelper)
      : m_Eulers(eulers), m_Phases(phases),
        m_CrystalStructures(crystalStructures), m_Mask(mask), m_Params(params),
        m_MasterValues(masterValues),
        m_MasterContributingCount(masterContributingCount),
        m_MasterTotalDeposits(masterTotalDeposits),
        m_MasterFailureCount(masterFailureCount), m_MergeMutex(mergeMutex),
        m_ShouldCancel(shouldCancel), m_ProgressHelper(progressHelper) {}

  void operator()(const Range &range) const {
    if (m_ShouldCancel) {
      return;
    }

    const std::size_t binCount = m_MasterValues.size();
    std::vector<double> localValues(binCount, 0.0);
    std::size_t localContributing = 0;
    std::size_t localTotalDeposits = 0;
    std::size_t localFailures = 0;

    auto progressMessenger = m_ProgressHelper.createProgressMessenger(
        std::chrono::milliseconds(500));
    const std::size_t crystalStructTuples =
        m_CrystalStructures.getNumberOfTuples();

    // Per-thread cache of EbsdLib's LaueOps table. The vector is cheap to copy
    // (vector of shared_ptr) and avoids re-querying the global table on each
    // voxel. Indexed by ebsdlib::CrystalStructure::*** codes.
    const std::vector<ebsdlib::LaueOps::Pointer> allOps =
        ebsdlib::LaueOps::GetAllOrientationOps();

    // Reusable scratch for the per-voxel symmetric expansion. Sized to the
    // maximum LaueOps count (24 for cubic) on first use; reused thereafter.
    std::vector<std::array<double, 3>> symEulers;

    for (std::size_t i = range.min(); i < range.max(); ++i) {
      // Cancel check kept at outer voxel loop only — inner per-symmetry-variant
      // work is short.
      if ((i & 0x3FFu) == 0 && m_ShouldCancel) {
        return;
      }

      // Mask filter (if present) and zero-phase filter.
      if (m_Mask != nullptr && !m_Mask->isTrue(i)) {
        progressMessenger.sendProgressMessage(1);
        continue;
      }
      const int32 phase = m_Phases[i];
      if (phase <= 0 ||
          static_cast<std::size_t>(phase) >= crystalStructTuples) {
        progressMessenger.sendProgressMessage(1);
        continue;
      }

      const uint32 code = m_CrystalStructures[static_cast<std::size_t>(phase)];
      if (code >= allOps.size() || allOps[code] == nullptr) {
        ++localFailures;
        progressMessenger.sendProgressMessage(1);
        continue;
      }
      const ebsdlib::LaueOps::Pointer laueOps = allOps[code];
      const std::size_t numOps = laueOps->getNumSymOps();

      // Promote single-precision EBSD inputs to double immediately.
      const ebsdlib::EulerDType inputEu(
          static_cast<double>(m_Eulers[3 * i + 0]),
          static_cast<double>(m_Eulers[3 * i + 1]),
          static_cast<double>(m_Eulers[3 * i + 2]));
      const ebsdlib::Matrix3X3D Gpassive =
          inputEu.toOrientationMatrix().toGMatrix();

      symEulers.resize(numOps);

      std::size_t deposits = 0;
      try {
        for (std::size_t k = 0; k < numOps; ++k) {
          // EbsdLib stores the ACTIVE symmetry operator; transpose for passive.
          // Symmetric variant: R = O_passive * G_passive_input.
          const ebsdlib::Matrix3X3D Op_active = laueOps->getMatSymOpD(k);
          const ebsdlib::Matrix3X3D R = Op_active.transpose() * Gpassive;

          // Extract Euler back via canonical om2eu (handles degenerate-PHI
          // and the eps-snap correctly). Wrap into [0, 2π) / [0, π] / [0, 2π).
          const ebsdlib::OrientationMatrixDType om(R[0], R[1], R[2], R[3], R[4],
                                                   R[5], R[6], R[7], R[8]);
          const ebsdlib::EulerDType outEu = om.toEuler();

          // Pragmatic conversion: ODFBuilder::accumulate currently consumes a
          // vector of std::array<double,3>. Changing that signature is out of
          // scope for this refactor (it's a pure binning concern, separate
          // from orientation math).
          symEulers[k] = {outEu[0], outEu[1], outEu[2]};
        }

        mtrsim::accumulate(symEulers, m_Params, localValues);
        deposits = numOps;
      } catch (const std::exception &) {
        // Accumulator-size mismatch (the only ODFBuilder::accumulate failure
        // mode). Counted but skipped — the master merger surfaces a warning if
        // any failures occurred.
        ++localFailures;
        progressMessenger.sendProgressMessage(1);
        continue;
      }
      ++localContributing;
      localTotalDeposits += deposits;
      progressMessenger.sendProgressMessage(1);
    }

    // Merge into master under a mutex.
    {
      std::lock_guard<std::mutex> lock(m_MergeMutex);
      for (std::size_t b = 0; b < binCount; ++b) {
        m_MasterValues[b] += localValues[b];
      }
      m_MasterContributingCount += localContributing;
      m_MasterTotalDeposits += localTotalDeposits;
      m_MasterFailureCount += localFailures;
    }
  }

private:
  const Float32Array &m_Eulers;
  const Int32Array &m_Phases;
  const UInt32Array &m_CrystalStructures;
  const MaskCompareUtilities::MaskCompare *m_Mask;
  const mtrsim::ODFBuildParams &m_Params;
  std::vector<double> &m_MasterValues;
  std::size_t &m_MasterContributingCount;
  std::size_t &m_MasterTotalDeposits;
  std::size_t &m_MasterFailureCount;
  std::mutex &m_MergeMutex;
  const std::atomic_bool &m_ShouldCancel;
  ProgressMessageHelper &m_ProgressHelper;
};

} // namespace

// -----------------------------------------------------------------------------
ComputeODF::ComputeODF(DataStructure &dataStructure,
                       const IFilter::MessageHandler &mesgHandler,
                       const std::atomic_bool &shouldCancel,
                       ComputeODFInputValues *inputValues)
    : m_DataStructure(dataStructure), m_InputValues(inputValues),
      m_ShouldCancel(shouldCancel), m_MessageHandler(mesgHandler) {}

// -----------------------------------------------------------------------------
ComputeODF::~ComputeODF() noexcept = default;

// -----------------------------------------------------------------------------
Result<> ComputeODF::operator()() {
  const auto &eulerAngles = m_DataStructure.getDataRefAs<Float32Array>(
      m_InputValues->eulerAnglesPath);
  const auto &phases =
      m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->phasesPath);
  const auto &crystalStructures = m_DataStructure.getDataRefAs<UInt32Array>(
      m_InputValues->crystalStructuresPath);
  // We get the pointer to the Array instead of a reference because it might not
  // have been set because the bool "use_mask" might have been false, but we do
  // NOT want to try to get the array 'on demand' in the loop. That is a BAD
  // idea as is it really slow to do that. (10x slower).
  std::unique_ptr<MaskCompareUtilities::MaskCompare> maskArray = nullptr;
  if (m_InputValues->useMask) {
    try {
      maskArray = MaskCompareUtilities::InstantiateMaskCompare(
          m_DataStructure, m_InputValues->maskPath);
    } catch (const std::out_of_range &exception) {
      // This really should NOT be happening as the path was verified during
      // preflight BUT we may be calling this from somewhere else that is NOT
      // going through the normal nx::core::IFilter API of Preflight and Execute
      std::string message =
          fmt::format("Mask Array DataPath does not exist or is not of the "
                      "correct type (Bool | UInt8) {}",
                      m_InputValues->maskPath.toString());
      return MakeErrorResult(-506, message);
    }
  }

  const DataPath outArrayPath =
      m_InputValues->outputImageGeometry
          .createChildPath(m_InputValues->cellAttrMatName)
          .createChildPath(m_InputValues->componentName);
  auto &outArray = m_DataStructure.getDataRefAs<Float64Array>(outArrayPath);
  auto &outStore = outArray.getDataStoreRef();

  const std::size_t binCount = static_cast<std::size_t>(m_InputValues->nphi1) *
                               static_cast<std::size_t>(m_InputValues->nPHI) *
                               static_cast<std::size_t>(m_InputValues->nphi2);
  if (outStore.getSize() != binCount) {
    return MakeErrorResult(
        -12210, fmt::format("Output ODF array size mismatch: have {} but "
                            "algorithm expected {} (= nphi1 * nPHI * nphi2).",
                            outStore.getSize(), binCount));
  }

  const std::size_t numVoxels = eulerAngles.getNumberOfTuples();

  m_MessageHandler(
      IFilter::Message::Type::Info,
      fmt::format("Computing ODF over {} voxels into {} bins ({} x {} x {}); "
                  "smoothing = {}.",
                  numVoxels, binCount, m_InputValues->nphi1,
                  m_InputValues->nPHI, m_InputValues->nphi2,
                  m_InputValues->applySmoothing ? "enabled" : "disabled"));

  // Set up the message helper / progress messenger BEFORE the parallel section
  // so the worker can stamp throttled per-voxel progress through its private
  // ProgressMessenger.
  MessageHelper messageHelper(m_MessageHandler);
  auto progressHelper = messageHelper.createProgressMessageHelper();
  progressHelper.setMaxProgresss(numVoxels);
  progressHelper.setProgressMessageTemplate("Compute ODF: {:.1f}%");

  std::vector<double> masterValues(binCount, 0.0);
  std::size_t masterContributingCount = 0;
  std::size_t masterTotalDeposits = 0;
  std::size_t masterFailureCount = 0;
  std::mutex mergeMutex;

  const mtrsim::ODFBuildParams params{
      m_InputValues->nphi1, m_InputValues->nPHI, m_InputValues->nphi2,
      m_InputValues->binSizeDeg, m_InputValues->applySmoothing};

  if (numVoxels == 0) {
    m_MessageHandler(IFilter::Message::Type::Info,
                     "No input voxels — output ODF array left at all zeros.");
    std::copy(masterValues.begin(), masterValues.end(), outStore.begin());
    return {};
  }

  ParallelDataAlgorithm parallelAlgorithm;
  parallelAlgorithm.setRange(0, numVoxels);
  parallelAlgorithm.execute(AccumulateWorker(
      eulerAngles, phases, crystalStructures, maskArray.get(), params,
      masterValues, masterContributingCount, masterTotalDeposits,
      masterFailureCount, mergeMutex, m_ShouldCancel, progressHelper));

  if (m_ShouldCancel) {
    return {};
  }

  // Normalize by the total count of symmetric-equivalent orientation deposits
  // (i.e. for each contributing voxel, the number of symmetric variants its
  // phase produces). This matches the MATLAB calc_ODF.m convention (line 80: N
  // = size(phi1_vec, 1) where phi1_vec is the post-symmetry-expansion
  // orientation list). The contributing-voxel count is retained as a semantic
  // guard for the "no voxels contributed" early-return path.
  if (masterContributingCount == 0) {
    m_MessageHandler(
        IFilter::Message::Type::Warning,
        "No voxels contributed to the ODF (mask filtered everything out, or "
        "all phases were 0). "
        "Output ODF array left at all zeros to avoid divide-by-zero.");
  } else {
    mtrsim::normalize(masterValues, static_cast<double>(masterTotalDeposits));
  }

  // Output-units conversion. The accumulator above holds raw count-density:
  // a histogram normalized so all Bunge bins sum to 1.0. That treats every bin
  // as if it has equal SO(3) volume, which is wrong — the SO(3) volume element
  // is sin(PHI) dphi1 dPHI dphi2, so bins near the PHI=0/180 poles cover
  // dramatically less orientation space than equatorial bins.
  //
  // For MUD output, divide each bin by its SO(3) bin volume and the uniform
  // SO(3) density 1/(8*pi^2):
  //   MUD[i,j,k] = acc[i,j,k] * 8*pi^2 / (step^3 * sin(PHI_center(j)))
  // This matches MTEX's plot(odf) convention: sub-uniform regions read < 1
  // MUD, peak texture reads several MUD.
  if (m_InputValues->outputUnits == 1ULL && masterContributingCount != 0) {
    const double stepRad = m_InputValues->binSizeDeg * std::numbers::pi / 180.0;
    const double scale = 8.0 * std::numbers::pi * std::numbers::pi /
                         (stepRad * stepRad * stepRad);
    for (int32_t j = 0; j < m_InputValues->nPHI; ++j) {
      const double phiCenter = (static_cast<double>(j) + 0.5) * stepRad;
      const double rowMul = scale / std::sin(phiCenter);
      for (int32_t i = 0; i < m_InputValues->nphi1; ++i) {
        const std::size_t baseIdx =
            (static_cast<std::size_t>(i) *
                 static_cast<std::size_t>(m_InputValues->nPHI) +
             static_cast<std::size_t>(j)) *
            static_cast<std::size_t>(m_InputValues->nphi2);
        for (int32_t k = 0; k < m_InputValues->nphi2; ++k) {
          masterValues[baseIdx + static_cast<std::size_t>(k)] *= rowMul;
        }
      }
    }
    m_MessageHandler(
        IFilter::Message::Type::Info,
        "Output units: MUD (per-PHI-row sin(PHI) Jacobian applied).");
  } else {
    m_MessageHandler(IFilter::Message::Type::Info,
                     "Output units: Count-Density (raw normalized histogram).");
  }

  m_MessageHandler(IFilter::Message::Type::Info,
                   fmt::format("ODF accumulation complete: {} voxel(s) "
                               "contributed; copying values to output array.",
                               masterContributingCount));

  // Bulk copy into the output Float64 store via iterators (matches the T2
  // polish convention).
  std::copy(masterValues.begin(), masterValues.end(), outStore.begin());

  // Surface per-voxel expansion failures (unknown crystal-structure code,
  // accumulator size mismatch) as a warning on the result so pipeline UIs can
  // show them. The run itself succeeds — the failing voxels are simply excluded
  // from the ODF.
  if (masterFailureCount > 0) {
    return MakeWarningVoidResult(
        -12213, fmt::format("Skipped {} voxel(s) due to expansion errors (e.g. "
                            "unknown crystal-structure code); those voxels did "
                            "not contribute to the ODF.",
                            masterFailureCount));
  }

  return {};
}

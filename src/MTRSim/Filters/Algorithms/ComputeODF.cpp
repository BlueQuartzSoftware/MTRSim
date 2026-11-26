#include "ComputeODF.hpp"

#include "simplnx/Common/Range.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/Utilities/MaskCompareUtilities.hpp"
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
#include <cmath>
#include <cstddef>
#include <exception>
#include <mutex>
#include <numbers>
#include <vector>

using namespace nx::core;

namespace
{

// Each worker invocation owns a private accumulator and contributing-voxel count.
// At completion, the worker merges its results into the master state under a mutex.
// This prevents concurrent writes to the master accumulator and counts.
//
// All orientation math goes through EbsdLib's `Euler<double>`,
// `OrientationMatrix<double>`, and `LaueOps` directly.
// The worker promotes Float32 EBSD inputs to double before orientation math.
// Double precision matches EbsdLib's convention and the project's orientation rule.
class AccumulateWorker
{
public:
  AccumulateWorker(const Float32Array& eulers, const Int32Array& phases, const UInt32Array& crystalStructures, int32 selectedPhase, const MaskCompareUtilities::MaskCompare* mask,
                   const mtrsim::ODFBuildParams& params, std::vector<double>& masterValues, std::size_t& masterContributingCount, std::size_t& masterTotalDeposits, std::size_t& masterFailureCount,
                   std::mutex& mergeMutex, const std::atomic_bool& shouldCancel, ComputeODF* algorithm)
  : m_Eulers(eulers)
  , m_Phases(phases)
  , m_CrystalStructures(crystalStructures)
  , m_SelectedPhase(selectedPhase)
  , m_Mask(mask)
  , m_Params(params)
  , m_MasterValues(masterValues)
  , m_MasterContributingCount(masterContributingCount)
  , m_MasterTotalDeposits(masterTotalDeposits)
  , m_MasterFailureCount(masterFailureCount)
  , m_MergeMutex(mergeMutex)
  , m_ShouldCancel(shouldCancel)
  , m_Algorithm(algorithm)
  {
  }

  void operator()(const Range& range) const
  {
    if(m_ShouldCancel)
    {
      return;
    }

    const std::size_t binCount = m_MasterValues.size();
    std::vector<double> localValues(binCount, 0.0);
    std::size_t localContributing = 0;
    std::size_t localTotalDeposits = 0;
    std::size_t localFailures = 0;

    const std::size_t progressIncrement = std::max((range.max() - range.min()) / 100, std::size_t{1});
    std::size_t pendingProgress = 0;
    const auto reportProgress = [&]() {
      pendingProgress++;
      if(pendingProgress >= progressIncrement)
      {
        m_Algorithm->sendThreadSafeProgressMessage(pendingProgress);
        pendingProgress = 0;
      }
    };
    const std::size_t crystalStructTuples = m_CrystalStructures.getNumberOfTuples();

    // Per-thread cache of EbsdLib's LaueOps table. The vector is cheap to copy
    // (vector of shared_ptr) and avoids re-querying the global table on each
    // voxel. Indexed by ebsdlib::CrystalStructure::*** codes.
    const std::vector<ebsdlib::LaueOps::Pointer> allOps = ebsdlib::LaueOps::GetAllOrientationOps();

    // Reusable scratch for the per-voxel symmetric expansion. Sized to the
    // maximum LaueOps count (24 for cubic) on first use; reused thereafter.
    std::vector<std::array<double, 3>> symEulers;

    for(std::size_t i = range.min(); i < range.max(); ++i)
    {
      // Cancel check kept at outer voxel loop only — inner per-symmetry-variant
      // work is short.
      if((i & 0x3FFu) == 0 && m_ShouldCancel)
      {
        return;
      }

      const int32 phase = m_Phases[i];
      if(phase != m_SelectedPhase || (m_Mask != nullptr && !m_Mask->isTrue(i)))
      {
        reportProgress();
        continue;
      }
      if(phase <= 0 || static_cast<std::size_t>(phase) >= crystalStructTuples)
      {
        reportProgress();
        continue;
      }

      const uint32 code = m_CrystalStructures[static_cast<std::size_t>(phase)];
      if(code >= allOps.size() || allOps[code] == nullptr)
      {
        ++localFailures;
        reportProgress();
        continue;
      }
      const ebsdlib::LaueOps::Pointer& laueOps = allOps[code];
      const std::size_t numOps = laueOps->getNumSymOps();

      // Promote single-precision EBSD inputs to double immediately.
      const ebsdlib::EulerDType inputEu(static_cast<double>(m_Eulers[3 * i + 0]), static_cast<double>(m_Eulers[3 * i + 1]), static_cast<double>(m_Eulers[3 * i + 2]));
      const ebsdlib::Matrix3X3D Gpassive = inputEu.toOrientationMatrix().toGMatrix();

      symEulers.resize(numOps);

      std::size_t deposits = 0;
      try
      {
        for(std::size_t k = 0; k < numOps; ++k)
        {
          // EbsdLib stores the ACTIVE symmetry operator; transpose for passive.
          // Symmetric variant: R = O_passive * G_passive_input.
          const ebsdlib::Matrix3X3D Op_active = laueOps->getMatSymOpD(k);
          const ebsdlib::Matrix3X3D R = Op_active.transpose() * Gpassive;

          // Extract Euler back via canonical om2eu (handles degenerate-PHI
          // and the eps-snap correctly). Wrap into [0, 2π) / [0, π] / [0, 2π).
          const ebsdlib::OrientationMatrixDType om(R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
          const ebsdlib::EulerDType outEu = om.toEuler();

          // Pragmatic conversion: ODFBuilder::accumulate currently consumes a
          // vector of std::array<double,3>. Changing that signature is out of
          // scope for this refactor (it's a pure binning concern, separate
          // from orientation math).
          symEulers[k] = {outEu[0], outEu[1], outEu[2]};
        }

        mtrsim::accumulate(symEulers, m_Params, localValues);
        deposits = numOps;
      } catch(const std::exception&)
      {
        // Accumulator-size mismatch (the only ODFBuilder::accumulate failure
        // mode). Counted but skipped — the master merger surfaces a warning if
        // any failures occurred.
        ++localFailures;
        reportProgress();
        continue;
      }
      ++localContributing;
      localTotalDeposits += deposits;
      reportProgress();
    }
    if(pendingProgress > 0)
    {
      m_Algorithm->sendThreadSafeProgressMessage(pendingProgress);
    }

    // Merge into master under a mutex.
    {
      std::lock_guard<std::mutex> lock(m_MergeMutex);
      for(std::size_t b = 0; b < binCount; ++b)
      {
        m_MasterValues[b] += localValues[b];
      }
      m_MasterContributingCount += localContributing;
      m_MasterTotalDeposits += localTotalDeposits;
      m_MasterFailureCount += localFailures;
    }
  }

private:
  const Float32Array& m_Eulers;
  const Int32Array& m_Phases;
  const UInt32Array& m_CrystalStructures;
  const int32 m_SelectedPhase;
  const MaskCompareUtilities::MaskCompare* m_Mask;
  const mtrsim::ODFBuildParams& m_Params;
  std::vector<double>& m_MasterValues;
  std::size_t& m_MasterContributingCount;
  std::size_t& m_MasterTotalDeposits;
  std::size_t& m_MasterFailureCount;
  std::mutex& m_MergeMutex;
  const std::atomic_bool& m_ShouldCancel;
  ComputeODF* m_Algorithm = nullptr;
};

} // namespace

// -----------------------------------------------------------------------------
ComputeODF::ComputeODF(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, ComputeODFInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
, m_ProgressThrottle(m_MessageHandler)
{
}

// -----------------------------------------------------------------------------
ComputeODF::~ComputeODF() noexcept = default;

void ComputeODF::sendThreadSafeProgressMessage(usize completedVoxels)
{
  const std::lock_guard<std::mutex> guard(m_ProgressMessageMutex);
  m_ProgressThrottle.incrementPercent(completedVoxels, 1);
}

// -----------------------------------------------------------------------------
Result<> ComputeODF::operator()()
{
  if(m_InputValues->contributingVoxelCount != nullptr)
  {
    *m_InputValues->contributingVoxelCount = 0;
  }
  const auto& eulerAngles = m_DataStructure.getDataRefAs<Float32Array>(m_InputValues->eulerAnglesPath);
  const auto& phases = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->phasesPath);
  const auto& crystalStructures = m_DataStructure.getDataRefAs<UInt32Array>(m_InputValues->crystalStructuresPath);
  const int32 selectedPhase = m_InputValues->selectedPhase;
  if(selectedPhase <= 0 || static_cast<usize>(selectedPhase) >= crystalStructures.getNumberOfTuples())
  {
    return MakeErrorResult(-12216, fmt::format("Phase Index ({}) must be positive and less than the tuple count ({}) of Crystal Structures '{}'.", selectedPhase, crystalStructures.getNumberOfTuples(),
                                               m_InputValues->crystalStructuresPath.toString()));
  }
  if(crystalStructures.getIDataStoreRef().getStoreType() == IDataStore::StoreType::Empty)
  {
    return MakeErrorResult(-12217, fmt::format("Crystal Structures '{}' values are not available for phase {}. Execute the input filter before computing the ODF.",
                                               m_InputValues->crystalStructuresPath.toString(), selectedPhase));
  }
  const uint32 crystalCode = crystalStructures[static_cast<usize>(selectedPhase)];
  const auto allOps = ebsdlib::LaueOps::GetAllOrientationOps();
  if(crystalCode >= allOps.size() || allOps[crystalCode] == nullptr)
  {
    return MakeErrorResult(-12217, fmt::format("Crystal Structures '{}' contains unsupported crystal-structure value {} at phase {}. Select a code from 0 through {}.",
                                               m_InputValues->crystalStructuresPath.toString(), crystalCode, selectedPhase, allOps.size() - 1));
  }
  // Resolve the optional mask once because repeated lookups in the voxel loop are expensive.
  std::unique_ptr<MaskCompareUtilities::MaskCompare> maskArray = nullptr;
  if(m_InputValues->useMask)
  {
    try
    {
      maskArray = MaskCompareUtilities::InstantiateMaskCompare(m_DataStructure, m_InputValues->maskPath);
    } catch(const std::out_of_range& exception)
    {
      // Direct algorithm callers can bypass the filter's parameter validation.
      std::string message = fmt::format("Mask Array DataPath does not exist or is not of the "
                                        "correct type (Bool | UInt8) {}",
                                        m_InputValues->maskPath.toString());
      return MakeErrorResult(-506, message);
    }
  }

  const DataPath outArrayPath = m_InputValues->outputImageGeometry.createChildPath(m_InputValues->cellAttrMatName).createChildPath(m_InputValues->componentName);
  auto& outArray = m_DataStructure.getDataRefAs<Float64Array>(outArrayPath);
  auto& outStore = outArray.getDataStoreRef();

  const std::size_t binCount = static_cast<std::size_t>(m_InputValues->nphi1) * static_cast<std::size_t>(m_InputValues->nPHI) * static_cast<std::size_t>(m_InputValues->nphi2);
  if(outStore.getSize() != binCount)
  {
    return MakeErrorResult(-12210, fmt::format("Output ODF array size mismatch: have {} but "
                                               "algorithm expected {} (= nphi1 * nPHI * nphi2).",
                                               outStore.getSize(), binCount));
  }

  const std::size_t numVoxels = eulerAngles.getNumberOfTuples();

  m_MessageHandler.sendInfoMessage(fmt::format("Computing ODF for phase {} over {} voxels into {} bins ({} x {} x {}); "
                                               "smoothing = {}.",
                                               selectedPhase, numVoxels, binCount, m_InputValues->nphi1, m_InputValues->nPHI, m_InputValues->nphi2,
                                               m_InputValues->applySmoothing ? "enabled" : "disabled"));

  m_ProgressThrottle.reset(numVoxels, fmt::format("Compute ODF for phase {}", selectedPhase));

  std::vector<double> masterValues(binCount, 0.0);
  std::size_t masterContributingCount = 0;
  std::size_t masterTotalDeposits = 0;
  std::size_t masterFailureCount = 0;
  std::mutex mergeMutex;

  const mtrsim::ODFBuildParams params{
      .nphi1 = m_InputValues->nphi1, .nPHI = m_InputValues->nPHI, .nphi2 = m_InputValues->nphi2, .binSizeDeg = m_InputValues->binSizeDeg, .smoothing = m_InputValues->applySmoothing};

  if(numVoxels == 0)
  {
    m_MessageHandler.sendInfoMessage("No input voxels — output ODF array left at all zeros.");
    std::copy(masterValues.begin(), masterValues.end(), outStore.begin());
    return {};
  }

  ParallelDataAlgorithm parallelAlgorithm;
  parallelAlgorithm.setRange(0, numVoxels);
  parallelAlgorithm.execute(AccumulateWorker(eulerAngles, phases, crystalStructures, selectedPhase, maskArray.get(), params, masterValues, masterContributingCount, masterTotalDeposits,
                                             masterFailureCount, mergeMutex, m_ShouldCancel, this));

  if(m_ShouldCancel)
  {
    return {};
  }

  if(m_InputValues->contributingVoxelCount != nullptr)
  {
    *m_InputValues->contributingVoxelCount = masterContributingCount;
  }

  // Normalize by the total count of symmetric-equivalent orientation deposits
  // (i.e. for each contributing voxel, the number of symmetric variants its
  // phase produces). This matches the MATLAB calc_ODF.m convention (line 80: N
  // = size(phi1_vec, 1) where phi1_vec is the post-symmetry-expansion
  // orientation list). The contributing-voxel count is retained as a semantic
  // guard for the "no voxels contributed" early-return path.
  if(masterContributingCount == 0)
  {
    m_MessageHandler.sendWarningMessage(fmt::format("No unmasked voxels of phase {} contributed to the ODF. Output ODF values are all zero.", selectedPhase));
  }
  else
  {
    mtrsim::normalize(masterValues, static_cast<double>(masterTotalDeposits));
  }

  // Count-density sums to 1 over the Bunge bins. The SO(3) volume element is sin(PHI) dphi1 dPHI dphi2.
  // Bins near the PHI poles cover less orientation space than equatorial bins.
  //
  // MUD divides count-density by the bin volume and the uniform SO(3) density, 1/(8*pi^2).
  //   MUD[i,j,k] = acc[i,j,k] * 8*pi^2 / (step^3 * sin(PHI_center(j)))
  // This matches MTEX's plot(odf) convention. Sub-uniform regions read below 1 MUD; texture peaks read above 1 MUD.
  if(m_InputValues->outputUnits == 1ULL && masterContributingCount != 0)
  {
    const double stepRad = m_InputValues->binSizeDeg * std::numbers::pi / 180.0;
    const double scale = 8.0 * std::numbers::pi * std::numbers::pi / (stepRad * stepRad * stepRad);
    for(int32_t j = 0; j < m_InputValues->nPHI; ++j)
    {
      const double phiCenter = (static_cast<double>(j) + 0.5) * stepRad;
      const double rowMul = scale / std::sin(phiCenter);
      for(int32_t i = 0; i < m_InputValues->nphi1; ++i)
      {
        const std::size_t baseIdx = (static_cast<std::size_t>(i) * static_cast<std::size_t>(m_InputValues->nPHI) + static_cast<std::size_t>(j)) * static_cast<std::size_t>(m_InputValues->nphi2);
        for(int32_t k = 0; k < m_InputValues->nphi2; ++k)
        {
          masterValues[baseIdx + static_cast<std::size_t>(k)] *= rowMul;
        }
      }
    }
    m_MessageHandler.sendInfoMessage("Output units: MUD (per-PHI-row sin(PHI) Jacobian applied).");
  }
  else
  {
    m_MessageHandler.sendInfoMessage("Output units: Count-Density (raw normalized histogram).");
  }

  m_MessageHandler.sendInfoMessage(fmt::format("ODF accumulation for phase {} complete: {} voxel(s) "
                                               "contributed; copying values to output array.",
                                               selectedPhase, masterContributingCount));

  // Bulk copy into the output Float64 store via iterators (matches the T2
  // polish convention).
  std::copy(masterValues.begin(), masterValues.end(), outStore.begin());

  // Surface per-voxel expansion failures (unknown crystal-structure code,
  // accumulator size mismatch) as a warning on the result so pipeline UIs can
  // show them. The run itself succeeds — the failing voxels are simply excluded
  // from the ODF.
  if(masterFailureCount > 0)
  {
    return MakeWarningVoidResult(-12213, fmt::format("Skipped {} voxel(s) due to expansion errors (e.g. "
                                                     "unknown crystal-structure code); those voxels did "
                                                     "not contribute to the ODF.",
                                                     masterFailureCount));
  }

  return {};
}

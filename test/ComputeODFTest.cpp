/**
 * @file
 * @brief Tests selected-phase ODF computation, metadata persistence, and MATLAB reference agreement.
 */

#include <algorithm>
#include <catch2/catch.hpp>

#include "MTRSim/Filters/Algorithms/ComputeODF.hpp"
#include "MTRSim/Filters/ComputeODFFilter.hpp"
#include "MTRSim/MTRSim_test_dirs.hpp"

#include "LibMTRSim/ODFFileIO.hpp"

#include "SimplnxCore/Filters/ReadDREAM3DFilter.hpp"
#include "SimplnxCore/Filters/WriteDREAM3DFilter.hpp"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataGroup.hpp"
#include "simplnx/DataStructure/DataStore.hpp"
#include "simplnx/DataStructure/EmptyStringStore.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/StringArray.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/Dream3dImportParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Pipeline/Pipeline.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <fmt/format.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <mutex>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

using namespace nx::core;
using namespace nx::core::unit_test;

namespace
{
// EbsdLib crystal-structure codes (matches ebsdlib::CrystalStructure
// namespace).
constexpr uint32_t k_HexagonalHigh = 0;
constexpr uint32_t k_CubicHigh = 1;

// Output paths used across tests.
const DataPath k_OutputGeomPath({"ODF"});
const std::string k_CellAttrMatName = "Cell Data";
const std::string k_ComponentName = "Component 1";

// Convenience: returns the linear index for the Float64 ODF cell-data array
// given (i_phi1, i_PHI, i_phi2). Mirrors the row-major linearization used
// inside ODFBuilder (phi1 slowest, phi2 fastest).
[[maybe_unused]] usize OdfLinear(int32_t iPhi1, int32_t iPHI, int32_t iPhi2, int32_t nPHI, int32_t nphi2)
{
  return static_cast<usize>(iPhi1) * static_cast<usize>(nPHI) * static_cast<usize>(nphi2) + static_cast<usize>(iPHI) * static_cast<usize>(nphi2) + static_cast<usize>(iPhi2);
}

// The fixture provides cell inputs and matching ensemble crystal codes and material names.
struct EbsdInputs
{
  DataStructure data_structure;
  DataPath ebsd_geom_path;
  DataPath cell_attr_mat_path;
  DataPath ensemble_attr_mat_path;
  DataPath euler_angles_path;
  DataPath phases_path;
  DataPath crystal_structures_path;
  DataPath material_names_path;
  DataPath mask_path;
};

EbsdInputs BuildSyntheticEbsd(const std::vector<std::array<float32, 3>>& eulersDeg, const std::vector<int32>& phaseLabels, uint32_t crystalCodeForPhase1, const std::vector<bool>* maskValues = nullptr,
                              bool includeCubicPhase2 = false)
{
  REQUIRE(eulersDeg.size() == phaseLabels.size());
  if(maskValues != nullptr)
  {
    REQUIRE(maskValues->size() == eulersDeg.size());
  }

  EbsdInputs out;
  const usize numVoxels = eulersDeg.size();

  // Lay out the input EBSD ImageGeom as a 1-D strip (numVoxels x 1 x 1). The
  // geometry shape is irrelevant to ODF output; only the per-voxel content
  // matters. Using {N, 1, 1} keeps the cell tuple-shape unambiguous so the
  // AttributeMatrix sees N tuples.
  ImageGeom* ebsdGeomPtr = ImageGeom::Create(out.data_structure, "EBSD");
  ebsdGeomPtr->setDimensions({numVoxels, 1, 1});
  out.ebsd_geom_path = DataPath({"EBSD"});

  AttributeMatrix* cellAmPtr = AttributeMatrix::Create(out.data_structure, k_CellAttrMatName, {numVoxels}, ebsdGeomPtr->getId());
  out.cell_attr_mat_path = out.ebsd_geom_path.createChildPath(k_CellAttrMatName);

  Float32ArrayType* eulerArrPtr = Float32ArrayType::CreateWithStore<DataStore<float32>>(out.data_structure, "EulerAngles", {numVoxels}, {3}, cellAmPtr->getId());
  Int32ArrayType* phaseArrPtr = Int32ArrayType::CreateWithStore<DataStore<int32>>(out.data_structure, "Phases", {numVoxels}, {1}, cellAmPtr->getId());

  constexpr double k_DegToRad = std::numbers::pi / 180.0;
  for(usize i = 0; i < numVoxels; ++i)
  {
    (*eulerArrPtr)[3 * i + 0] = static_cast<float32>(eulersDeg[i][0] * k_DegToRad);
    (*eulerArrPtr)[3 * i + 1] = static_cast<float32>(eulersDeg[i][1] * k_DegToRad);
    (*eulerArrPtr)[3 * i + 2] = static_cast<float32>(eulersDeg[i][2] * k_DegToRad);
    (*phaseArrPtr)[i] = phaseLabels[i];
  }

  out.euler_angles_path = out.cell_attr_mat_path.createChildPath("EulerAngles");
  out.phases_path = out.cell_attr_mat_path.createChildPath("Phases");

  if(maskValues != nullptr)
  {
    BoolArrayType* maskArrPtr = BoolArrayType::CreateWithStore<DataStore<bool>>(out.data_structure, "Mask", {numVoxels}, {1}, cellAmPtr->getId());
    for(usize i = 0; i < numVoxels; ++i)
    {
      (*maskArrPtr)[i] = (*maskValues)[i];
    }
    out.mask_path = out.cell_attr_mat_path.createChildPath("Mask");
  }

  // Ensemble arrays can reside outside the input geometry.
  const usize ensembleCount = includeCubicPhase2 ? 3 : 2;
  AttributeMatrix* ensembleAmPtr = AttributeMatrix::Create(out.data_structure, "EnsembleData", {ensembleCount});
  out.ensemble_attr_mat_path = DataPath({"EnsembleData"});

  UInt32ArrayType* csArrPtr = UInt32ArrayType::CreateWithStore<DataStore<uint32>>(out.data_structure, "CrystalStructures", {ensembleCount}, {1}, ensembleAmPtr->getId());
  (*csArrPtr)[0] = 999; // phase 0 = "no phase" — never indexed by the algorithm
  (*csArrPtr)[1] = crystalCodeForPhase1;
  out.crystal_structures_path = out.ensemble_attr_mat_path.createChildPath("CrystalStructures");
  std::vector<std::string> materialNames = {"Unknown", "Phase 1"};
  if(includeCubicPhase2)
  {
    (*csArrPtr)[2] = k_CubicHigh;
    materialNames.emplace_back("Cubic Phase 2");
  }
  REQUIRE(StringArray::CreateWithValues(out.data_structure, "Material Names", {ensembleCount}, materialNames, ensembleAmPtr->getId()) != nullptr);
  out.material_names_path = out.ensemble_attr_mat_path.createChildPath("Material Names");

  return out;
}

// Pre-fills a base Arguments object with all the required input paths. Tests
// then override individual values as needed.
Arguments MakeBaseArgs(const EbsdInputs& inputs, bool applySmoothing, float32 binSizeDeg)
{
  Arguments args;
  args.insertOrAssign(ComputeODFFilter::k_SaveAsImageGeometry_Key, false);
  args.insertOrAssign(ComputeODFFilter::k_ApplySmoothing_Key, std::make_any<bool>(applySmoothing));
  args.insertOrAssign(ComputeODFFilter::k_BinSizeDeg_Key, std::make_any<float32>(binSizeDeg));
  args.insertOrAssign(ComputeODFFilter::k_EulerAngles_Key, std::make_any<DataPath>(inputs.euler_angles_path));
  args.insertOrAssign(ComputeODFFilter::k_Phases_Key, std::make_any<DataPath>(inputs.phases_path));
  args.insertOrAssign(ComputeODFFilter::k_CrystalStructures_Key, std::make_any<DataPath>(inputs.crystal_structures_path));
  args.insertOrAssign(ComputeODFFilter::k_MaterialNames_Key, inputs.material_names_path);
  args.insertOrAssign(ComputeODFFilter::k_PhaseIndex_Key, int32{1});
  args.insertOrAssign(ComputeODFFilter::k_UseMask_Key, std::make_any<bool>(false));
  args.insertOrAssign(ComputeODFFilter::k_MaskArrayPath_Key, std::make_any<DataPath>(DataPath{}));
  args.insertOrAssign(ComputeODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(k_OutputGeomPath));
  args.insertOrAssign(ComputeODFFilter::k_CellAttrMatName_Key, std::make_any<DataObjectNameParameter::ValueType>(k_CellAttrMatName));
  args.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::make_any<DataObjectNameParameter::ValueType>(k_ComponentName));
  return args;
}

// Returns the code of the first error in a preflight result, or 0 if there are
// no errors.
int32 FirstErrorCode(const IFilter::PreflightResult& r)
{
  const auto& errors = r.output_actions.errors();
  return errors.empty() ? 0 : errors[0].code;
}

// Returns (sum, nonZeroBinCount) over the output ODF array.
std::pair<double, usize> SumAndNonZeroCount(const DataStructure& dataStructure)
{
  const DataPath odfPath = k_OutputGeomPath.createChildPath(k_CellAttrMatName).createChildPath(k_ComponentName);
  const auto& odfArr = dataStructure.getDataRefAs<Float64ArrayType>(odfPath);
  const auto& store = odfArr.getDataStoreRef();
  double sum = 0.0;
  usize nonZero = 0;
  for(usize i = 0; i < store.getSize(); ++i)
  {
    const double v = store[i];
    sum += v;
    if(v != 0.0)
    {
      ++nonZero;
    }
  }
  return {sum, nonZero};
}

} // namespace

TEST_CASE("MTRSim::ComputeODFFilter: Single HCP voxel, no smoothing, produces "
          "a normalized ODF",
          "[MTRSim][ComputeODFFilter]")
{
  // 1 voxel of HCP at (10, 20, 30) deg. With smoothing off, each of the 12
  // symmetric variants deposits exactly 1.0 into a (possibly shared) bin.
  // Normalization divides by the total deposit count = 12 (matches MATLAB
  // calc_ODF.m), so the resulting sum should be exactly 1.0.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto [sum, nonZero] = SumAndNonZeroCount(inputs.data_structure);
  REQUIRE(sum == Approx(1.0).margin(1.0e-9));
  // Some symmetric variants may map to the same bin, so >=1 and <=12 nonzero
  // bins are valid.
  REQUIRE(nonZero >= 1);
  REQUIRE(nonZero <= 12);
}

TEST_CASE("MTRSim::ComputeODFFilter: Single Cubic voxel, no smoothing, "
          "produces a normalized ODF",
          "[MTRSim][ComputeODFFilter]")
{
  // 1 voxel of Cubic_High at (10, 20, 30) deg. 24 symmetric variants × 1.0
  // each, normalized by the total deposit count = 24 (matches MATLAB
  // calc_ODF.m), so the resulting sum should be exactly 1.0.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_CubicHigh);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto [sum, nonZero] = SumAndNonZeroCount(inputs.data_structure);
  REQUIRE(sum == Approx(1.0).margin(1.0e-9));
  REQUIRE(nonZero >= 1);
  REQUIRE(nonZero <= 24);
}

TEST_CASE("MTRSim::ComputeODFFilter: Smoothing distributes per MATLAB weights", "[MTRSim][ComputeODFFilter]")
{
  // 1 HCP voxel near a bin edge (12.5, 12.5, 12.5) deg. 12.5/5.0 = 2.5 so this
  // value sits exactly on the boundary between bins 2 and 3 (NOT at a bin
  // center). With smoothing on, each of the 12 symmetric variants
  // distributes 1.0 across 27 bins (the MATLAB tri-linear weights sum to 1.0
  // regardless of where inside the cube the sample falls). Total = 12 x 1.0;
  // normalize by total deposit count N=12 -> final sum = 1.0.
  EbsdInputs inputs = BuildSyntheticEbsd({{{12.5f, 12.5f, 12.5f}}}, {1}, k_HexagonalHigh);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/true, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto [sum, nonZero] = SumAndNonZeroCount(inputs.data_structure);
  REQUIRE(sum == Approx(1.0).margin(1.0e-9));
  // Smoothing across 27 bins per variant means many more nonzero bins are
  // expected.
  REQUIRE(nonZero >= 12);
}

TEST_CASE("MTRSim::ComputeODFFilter: MUD output applies per-row sin(PHI) Jacobian", "[MTRSim][ComputeODFFilter]")
{
  // Run the filter twice on identical input (same Eulers, smoothing on, default
  // hex symmetry) — once in Count-Density mode and once in MUD mode — then
  // verify that for every non-zero bin the per-bin ratio matches the
  // per-PHI-row Jacobian
  //   factor = 8 * pi^2 / (step^3 * sin(PHI_center(j)))
  // This locks the conversion formula in place. Any future change (wrong row
  // index, missing pi^2, swapped step exponent, etc.) breaks this immediately.
  //
  // The input mixes a low-PHI sample (PHI = 5 deg, large Jacobian factor) and a
  // high-PHI sample (PHI = 87.5 deg, near-equator factor) so both extremes of
  // the row multiplier are exercised.
  std::vector<std::array<float32, 3>> eulers = {
      {10.0f, 20.0f, 30.0f},
      {40.0f, 60.0f, 80.0f},
      {15.0f, 5.0f, 12.5f},
      {200.0f, 87.5f, 200.0f},
  };
  std::vector<int32> phases = {1, 1, 1, 1};

  EbsdInputs inputsA = BuildSyntheticEbsd(eulers, phases, k_HexagonalHigh);
  Arguments argsA = MakeBaseArgs(inputsA, /*applySmoothing=*/true, /*binSizeDeg=*/5.0f);
  argsA.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  ComputeODFFilter filterA;
  auto preflightA = filterA.preflight(inputsA.data_structure, argsA);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightA.output_actions);
  auto executeA = filterA.execute(inputsA.data_structure, argsA);
  SIMPLNX_RESULT_REQUIRE_VALID(executeA.result);

  EbsdInputs inputsB = BuildSyntheticEbsd(eulers, phases, k_HexagonalHigh);
  Arguments argsB = MakeBaseArgs(inputsB, /*applySmoothing=*/true, /*binSizeDeg=*/5.0f);
  argsB.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  ComputeODFFilter filterB;
  auto preflightB = filterB.preflight(inputsB.data_structure, argsB);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightB.output_actions);
  auto executeB = filterB.execute(inputsB.data_structure, argsB);
  SIMPLNX_RESULT_REQUIRE_VALID(executeB.result);

  const DataPath odfPath = k_OutputGeomPath.createChildPath(k_CellAttrMatName).createChildPath(k_ComponentName);
  const auto& cdStore = inputsA.data_structure.getDataRefAs<Float64ArrayType>(odfPath).getDataStoreRef();
  const auto& mudStore = inputsB.data_structure.getDataRefAs<Float64ArrayType>(odfPath).getDataStoreRef();
  REQUIRE(cdStore.getSize() == mudStore.getSize());

  constexpr int32_t k_Nphi1 = 72;
  constexpr int32_t k_NPhi = 36;
  constexpr int32_t k_Nphi2 = 72;
  REQUIRE(cdStore.getSize() == static_cast<usize>(k_Nphi1) * static_cast<usize>(k_NPhi) * static_cast<usize>(k_Nphi2));

  const double stepRad = 5.0 * std::numbers::pi / 180.0;
  const double scale = 8.0 * std::numbers::pi * std::numbers::pi / (stepRad * stepRad * stepRad);

  std::size_t binsCompared = 0;
  std::size_t binsFailed = 0;
  double maxRelErr = 0.0;
  for(int32_t j = 0; j < k_NPhi; ++j)
  {
    const double phiCenter = (static_cast<double>(j) + 0.5) * stepRad;
    const double expectedRowMul = scale / std::sin(phiCenter);
    for(int32_t i = 0; i < k_Nphi1; ++i)
    {
      for(int32_t k = 0; k < k_Nphi2; ++k)
      {
        const usize idx = OdfLinear(i, j, k, k_NPhi, k_Nphi2);
        const double cd = cdStore[idx];
        const double mud = mudStore[idx];
        if(cd == 0.0)
        {
          REQUIRE(mud == 0.0);
          continue;
        }
        ++binsCompared;
        const double expected = cd * expectedRowMul;
        const double rel = std::abs(mud - expected) / std::max(std::abs(expected), 1.0e-300);
        maxRelErr = std::max(rel, maxRelErr);
        if(rel > 1.0e-12)
        {
          ++binsFailed;
        }
      }
    }
  }
  INFO(fmt::format("compared {} non-zero bins; max relative error = {:.3e}; failing = {}", binsCompared, maxRelErr, binsFailed));
  REQUIRE(binsCompared > 0);
  REQUIRE(binsFailed == 0);
}

TEST_CASE("MTRSim::ComputeODFFilter: MUD integrates to 8*pi^2 on SO(3)", "[MTRSim][ComputeODFFilter]")
{
  // Single-number sanity check: by construction the count-density ODF sums
  // to 1.0 (it's a probability distribution on Bunge bins). Converting to MUD
  // multiplies each bin by 8*pi^2 / (step^3 * sin(PHI_center)). If we then
  // weight each bin by its true SO(3) volume element (step^3 *
  // sin(PHI_center)), we should recover
  //   sum_bins (MUD[bin] * step^3 * sin(PHI_center(bin))) = 8*pi^2
  // i.e. the SO(3) volume. This is the closed-form integral identity for MUD.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}, {{50.0f, 70.0f, 110.0f}}}, {1, 1}, k_HexagonalHigh);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/true, /*binSizeDeg=*/5.0f);
  args.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const DataPath odfPath = k_OutputGeomPath.createChildPath(k_CellAttrMatName).createChildPath(k_ComponentName);
  const auto& mudStore = inputs.data_structure.getDataRefAs<Float64ArrayType>(odfPath).getDataStoreRef();

  constexpr int32_t k_Nphi1 = 72;
  constexpr int32_t k_NPhi = 36;
  constexpr int32_t k_Nphi2 = 72;
  const double stepRad = 5.0 * std::numbers::pi / 180.0;
  const double step3 = stepRad * stepRad * stepRad;

  double integral = 0.0;
  for(int32_t j = 0; j < k_NPhi; ++j)
  {
    const double rowVol = step3 * std::sin((static_cast<double>(j) + 0.5) * stepRad);
    double rowSum = 0.0;
    for(int32_t i = 0; i < k_Nphi1; ++i)
    {
      for(int32_t k = 0; k < k_Nphi2; ++k)
      {
        rowSum += mudStore[OdfLinear(i, j, k, k_NPhi, k_Nphi2)];
      }
    }
    integral += rowSum * rowVol;
  }

  // 8*pi^2 ≈ 78.9568. Tolerance of 1e-9 is comfortable since both sides are
  // pure double-precision arithmetic with no symmetry-cancellation traps.
  const double expected = 8.0 * std::numbers::pi * std::numbers::pi;
  INFO(fmt::format("integral = {:.10f}, expected = {:.10f}", integral, expected));
  REQUIRE(integral == Approx(expected).margin(1.0e-9));
}

TEST_CASE("MTRSim::ComputeODFFilter: Mask filters voxels", "[MTRSim][ComputeODFFilter]")
{
  // 4 HCP voxels with distinct Eulers, mask = {true, false, true, false}.
  // 2 contributing voxels × 12 variants = 24 deposits; normalized by total
  // deposit count N=24 (matches MATLAB calc_ODF.m) → sum = 1.0.
  std::vector<std::array<float32, 3>> eulers = {{{10.0f, 20.0f, 30.0f}}, {{40.0f, 50.0f, 60.0f}}, {{70.0f, 80.0f, 90.0f}}, {{15.0f, 25.0f, 35.0f}}};
  std::vector<int32> phases = {1, 1, 1, 1};
  std::vector<bool> mask = {true, false, true, false};
  EbsdInputs inputs = BuildSyntheticEbsd(eulers, phases, k_HexagonalHigh, &mask);

  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  args.insertOrAssign(ComputeODFFilter::k_UseMask_Key, std::make_any<bool>(true));
  args.insertOrAssign(ComputeODFFilter::k_MaskArrayPath_Key, std::make_any<DataPath>(inputs.mask_path));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto [sum, nonZero] = SumAndNonZeroCount(inputs.data_structure);
  REQUIRE(sum == Approx(1.0).margin(1.0e-9));
  REQUIRE(nonZero >= 1);
}

TEST_CASE("MTRSim::ComputeODFFilter: Progress messages drive the numeric progress field", "[MTRSim][ComputeODFFilter][Messages]")
{
  std::vector<std::array<float32, 3>> eulers(4096, {10.0F, 20.0F, 30.0F});
  std::vector<int32> phases(eulers.size(), 1);
  EbsdInputs inputs = BuildSyntheticEbsd(eulers, phases, k_HexagonalHigh);
  Arguments args = MakeBaseArgs(inputs, false, 5.0F);

  std::mutex messagesMutex;
  std::vector<IFilter::Message> messages;
  const IFilter::MessageHandler handler{[&](const IFilter::Message& message) {
    const std::lock_guard<std::mutex> lock(messagesMutex);
    messages.push_back(message);
  }};

  ComputeODFFilter filter;
  const auto executeResult = filter.execute(inputs.data_structure, args, nullptr, handler);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  std::vector<IFilter::Message> progressMessages;
  for(const auto& message : messages)
  {
    if(message.type == IFilter::Message::Type::Progress)
    {
      progressMessages.push_back(message);
    }
  }
  REQUIRE_FALSE(progressMessages.empty());
  REQUIRE(progressMessages.back().progress == 100);
  REQUIRE(progressMessages.back().message.find("Compute ODF for phase 1") != std::string::npos);
}

TEST_CASE("MTRSim::ComputeODFFilter: Phase 0 voxels are skipped", "[MTRSim][ComputeODFFilter]")
{
  // 4 HCP voxels with phase pattern {1, 0, 1, 0}, all same Euler. No mask.
  // 2 contributing × 12 variants = 24 deposits; normalize by total deposit
  // count N=24 (matches MATLAB calc_ODF.m) → sum = 1.0.
  std::vector<std::array<float32, 3>> eulers = {{{10.0f, 20.0f, 30.0f}}, {{10.0f, 20.0f, 30.0f}}, {{10.0f, 20.0f, 30.0f}}, {{10.0f, 20.0f, 30.0f}}};
  std::vector<int32> phases = {1, 0, 1, 0};
  EbsdInputs inputs = BuildSyntheticEbsd(eulers, phases, k_HexagonalHigh);

  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto [sum, nonZero] = SumAndNonZeroCount(inputs.data_structure);
  REQUIRE(sum == Approx(1.0).margin(1.0e-9));
  REQUIRE(nonZero >= 1);
}

TEST_CASE("MTRSim::ComputeODFFilter: Empty mask DataPath disables the mask", "[MTRSim][ComputeODFFilter]")
{
  // 4 HCP voxels, all phase=1, distinct Eulers. UseMask = false (mask path is
  // empty). 4 contributing × 12 variants = 48 deposits; normalize by total
  // deposit count N=48 (matches MATLAB calc_ODF.m) → sum = 1.0.
  std::vector<std::array<float32, 3>> eulers = {{{10.0f, 20.0f, 30.0f}}, {{40.0f, 50.0f, 60.0f}}, {{70.0f, 80.0f, 90.0f}}, {{15.0f, 25.0f, 35.0f}}};
  std::vector<int32> phases = {1, 1, 1, 1};
  EbsdInputs inputs = BuildSyntheticEbsd(eulers, phases, k_HexagonalHigh);

  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto [sum, nonZero] = SumAndNonZeroCount(inputs.data_structure);
  REQUIRE(sum == Approx(1.0).margin(1.0e-9));
  REQUIRE(nonZero >= 1);
}

TEST_CASE("MTRSim::ComputeODFFilter: Preflight rejects bad bin size", "[MTRSim][ComputeODFFilter]")
{
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/7.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  REQUIRE(preflightResult.output_actions.invalid());
  REQUIRE(FirstErrorCode(preflightResult) == -12200);
}

TEST_CASE("MTRSim::ComputeODFFilter: Preflight rejects euler_angles wrong "
          "component count",
          "[MTRSim][ComputeODFFilter]")
{
  // Build a normal DataStructure, then swap in a wrong-component-count Float32
  // array for euler_angles. The ArraySelectionParameter component-shape
  // constraint will cause the IFilter::preflight wrapper to flag it before
  // preflightImpl is even called.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);

  // Add an alternate single-component Float32 array on the same cell AM.
  AttributeMatrix& cellAm = inputs.data_structure.getDataRefAs<AttributeMatrix>(inputs.cell_attr_mat_path);
  Float32ArrayType::CreateWithStore<DataStore<float32>>(inputs.data_structure, "BadEulers", {1}, {1}, cellAm.getId());
  const DataPath badPath = inputs.cell_attr_mat_path.createChildPath("BadEulers");

  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  args.insertOrAssign(ComputeODFFilter::k_EulerAngles_Key, std::make_any<DataPath>(badPath));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  REQUIRE(preflightResult.output_actions.invalid());
  // The ArraySelectionParameter component-shape constraint is enforced by the
  // IFilter::preflight wrapper (code -208 =
  // FilterParameter::Constants::k_Validate_TupleShapeValue) BEFORE
  // preflightImpl runs, so our own -12201 code never gets a chance to fire.
  // Assert the wrapper code explicitly.
  REQUIRE(FirstErrorCode(preflightResult) == -208);
}

TEST_CASE("MTRSim::ComputeODFFilter: Preflight rejects mismatched attribute matrices", "[MTRSim][ComputeODFFilter]")
{
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);

  // Create a second cell AttributeMatrix on the same geometry with its own
  // Phases array, so euler_angles and phases live on different parents.
  ImageGeom& ebsdGeom = inputs.data_structure.getDataRefAs<ImageGeom>(inputs.ebsd_geom_path);
  AttributeMatrix* otherAmPtr = AttributeMatrix::Create(inputs.data_structure, "OtherCellData", {1}, ebsdGeom.getId());
  Int32ArrayType* otherPhasesPtr = Int32ArrayType::CreateWithStore<DataStore<int32>>(inputs.data_structure, "Phases", {1}, {1}, otherAmPtr->getId());
  (*otherPhasesPtr)[0] = 1;
  const DataPath otherPhasesPath = inputs.ebsd_geom_path.createChildPath("OtherCellData").createChildPath("Phases");

  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  args.insertOrAssign(ComputeODFFilter::k_Phases_Key, std::make_any<DataPath>(otherPhasesPath));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  REQUIRE(preflightResult.output_actions.invalid());
  REQUIRE(FirstErrorCode(preflightResult) == -12205);
}

TEST_CASE("MTRSim::ComputeODFFilter: Preflight populates updatedValues", "[MTRSim][ComputeODFFilter]")
{
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/true, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);

  REQUIRE(!preflightResult.output_values.empty());

  const bool foundBinSize = std::any_of(preflightResult.output_values.begin(), preflightResult.output_values.end(), [](const IFilter::PreflightValue& v) { return v.name == "Bin Size [deg]"; });
  REQUIRE(foundBinSize);

  const bool foundBins = std::any_of(preflightResult.output_values.begin(), preflightResult.output_values.end(), [](const IFilter::PreflightValue& v) { return v.name == "Bins (phi1 x PHI x phi2)"; });
  REQUIRE(foundBins);
}

// -----------------------------------------------------------------------------
// Append-mode tests (Milestone AJ, Task 7b). Build on T7a: first run Create New
// to populate an ODF geometry, then run Append against that same geometry.
// -----------------------------------------------------------------------------

namespace
{
// Builds a pre-existing ImageGeom with the given XYZ dimensions, XYZ spacing,
// and a cell AttributeMatrix named "Cell Data". Optionally pre-populates a
// "Component 1" Float64 cell-data array with a constant value so tests can
// verify it's left alone by an Append call.
struct ExistingOdfGeom
{
  DataPath geom_path;
  DataPath cell_attr_mat_path;
};

ExistingOdfGeom BuildExistingOdfGeom(DataStructure& dataStructure, const std::array<usize, 3>& dimsXYZ, const std::array<float32, 3>& spacingXYZ, bool seedComponent1 = false, double seedValue = 0.0)
{
  ExistingOdfGeom out;
  ImageGeom* geomPtr = ImageGeom::Create(dataStructure, "ODF");
  geomPtr->setDimensions({dimsXYZ[0], dimsXYZ[1], dimsXYZ[2]});
  REQUIRE(geomPtr->setSpacing({spacingXYZ[0], spacingXYZ[1], spacingXYZ[2]}).valid());
  REQUIRE(geomPtr->setOrigin({0.0f, 0.0f, 0.0f}).valid());
  out.geom_path = DataPath({"ODF"});

  // Tuple shape is ZYX for cell-data arrays — match the convention used by
  // Create New mode.
  const std::vector<usize> tupleShapeZYX = {dimsXYZ[2], dimsXYZ[1], dimsXYZ[0]};
  AttributeMatrix* cellAmPtr = AttributeMatrix::Create(dataStructure, k_CellAttrMatName, tupleShapeZYX, geomPtr->getId());
  // Register the AttributeMatrix with the ImageGeom as its cell data so
  // getCellDataPath() works.
  geomPtr->setCellData(cellAmPtr->getId());
  out.cell_attr_mat_path = out.geom_path.createChildPath(k_CellAttrMatName);

  if(seedComponent1)
  {
    Float64ArrayType* seedArrPtr = Float64ArrayType::CreateWithStore<DataStore<float64>>(dataStructure, k_ComponentName, tupleShapeZYX, {1}, cellAmPtr->getId());
    auto& store = seedArrPtr->getDataStoreRef();
    for(usize i = 0; i < store.getSize(); ++i)
    {
      store[i] = seedValue;
    }
  }
  return out;
}

// Returns (sum, nonZeroCount) over an arbitrary Float64 array path.
std::pair<double, usize> SumAndNonZeroOf(const DataStructure& dataStructure, const DataPath& path)
{
  const auto& arr = dataStructure.getDataRefAs<Float64ArrayType>(path);
  const auto& store = arr.getDataStoreRef();
  double sum = 0.0;
  usize nonZero = 0;
  for(usize i = 0; i < store.getSize(); ++i)
  {
    const double v = store[i];
    sum += v;
    if(v != 0.0)
    {
      ++nonZero;
    }
  }
  return {sum, nonZero};
}
} // namespace

TEST_CASE("MTRSim::ComputeODFFilter: Selected phase controls symmetry and masked contributions", "[MTRSim][ComputeODFFilter]")
{
  const int32 selectedPhase = GENERATE(1, 2);
  const bool useMask = GENERATE(false, true);
  const std::vector<bool> mask = {true, false, true};
  auto inputs = BuildSyntheticEbsd({{{11.0f, 22.0f, 33.0f}}, {{41.0f, 57.0f, 83.0f}}, {{17.0f, 43.0f, 79.0f}}}, {1, 1, 2}, k_HexagonalHigh, &mask, true);
  auto args = MakeBaseArgs(inputs, false, 5.0F);
  args.insertOrAssign("phase_index", selectedPhase);
  args.insertOrAssign(ComputeODFFilter::k_UseMask_Key, useMask);
  args.insertOrAssign(ComputeODFFilter::k_MaskArrayPath_Key, inputs.mask_path);
  const auto result = ComputeODFFilter().execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(result.result);
  const auto [sum, nonZero] = SumAndNonZeroCount(inputs.data_structure);
  REQUIRE(sum == Approx(1.0).margin(1.0e-12));
  // These general orientations have distinct 5-degree bins for each symmetry operator.
  REQUIRE(nonZero == (selectedPhase == 1 && useMask ? 12 : 24));
  const DataPath componentPath({"ODF", "Cell Data", "Component 1"});
  REQUIRE_NOTHROW(inputs.data_structure.getDataRefAs<Float64ArrayType>(componentPath));
  const auto& values = inputs.data_structure.getDataRefAs<Float64ArrayType>(componentPath);
  // The identity operator places the phase-1 orientation at (11, 22, 33) degrees.
  REQUIRE(values[OdfLinear(2, 4, 6, 36, 72)] == Approx(selectedPhase == 2 ? 0.0 : (useMask ? 1.0 / 12.0 : 1.0 / 24.0)).margin(1.0e-12));
}

TEST_CASE("MTRSim::ComputeODFFilter: Units metadata is initialized by preflight actions", "[MTRSim][ComputeODFFilter]")
{
  const int outputSetup = GENERATE(0, 1, 2);
  const ChoicesParameter::ValueType units = GENERATE(0ULL, 1ULL);
  const bool pipelinePreflight = GENERATE(false, true);
  auto inputs = BuildSyntheticEbsd({{{11.0F, 22.0F, 33.0F}}}, {1}, k_HexagonalHigh);
  auto args = MakeBaseArgs(inputs, false, 5.0F);
  args.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, units);
  if(outputSetup != 0)
  {
    BuildExistingOdfGeom(inputs.data_structure, {72, 36, 72}, {5.0F, 5.0F, 5.0F});
    if(outputSetup == 2)
    {
      const auto* geometryPtr = inputs.data_structure.getDataAs<ImageGeom>(k_OutputGeomPath);
      REQUIRE(geometryPtr != nullptr);
      REQUIRE(DataGroup::Create(inputs.data_structure, "ODF Metadata", geometryPtr->getId()) != nullptr);
    }
    args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, ChoicesParameter::ValueType{1});
    args.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, k_OutputGeomPath);
  }
  if(pipelinePreflight)
  {
    Pipeline pipeline;
    REQUIRE(pipeline.push_back(std::make_unique<ComputeODFFilter>(), args));
    const std::atomic_bool shouldCancel = false;
    REQUIRE(pipeline.preflight(inputs.data_structure, shouldCancel));
  }
  else
  {
    const auto preflight = ComputeODFFilter().preflight(inputs.data_structure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflight.output_actions);
    const auto applied = preflight.output_actions.value().applyAll(inputs.data_structure, IDataAction::Mode::Preflight);
    SIMPLNX_RESULT_REQUIRE_VALID(applied);
  }
  const auto* unitsArrayPtr = inputs.data_structure.getDataAs<StringArray>(DataPath({"ODF", "ODF Metadata", "Component 1 Units"}));
  REQUIRE(unitsArrayPtr != nullptr);
  REQUIRE(unitsArrayPtr->getNumberOfTuples() == 1);
  REQUIRE((*unitsArrayPtr)[0] == (units == 0 ? "Count-Density" : "MUD"));
}

TEST_CASE("MTRSim::ComputeODFFilter: Rejects invalid ensemble selections", "[MTRSim][ComputeODFFilter]")
{
  auto inputs = BuildSyntheticEbsd({{{11.0F, 22.0F, 33.0F}}}, {1}, k_HexagonalHigh);
  auto args = MakeBaseArgs(inputs, false, 5.0F);
  int32 expectedCode = -12216;
  SECTION("zero phase")
  {
    args.insertOrAssign("phase_index", int32{0});
  }
  SECTION("negative phase")
  {
    args.insertOrAssign("phase_index", int32{-1});
  }
  SECTION("phase outside ensemble")
  {
    args.insertOrAssign("phase_index", int32{2});
  }
  SECTION("wrong material object type")
  {
    args.insertOrAssign("material_names_path", inputs.crystal_structures_path);
    expectedCode = -12214;
  }
  SECTION("different ensemble tuple counts")
  {
    REQUIRE(StringArray::CreateWithValues(inputs.data_structure, "Short Names", {1}, {"Only"}) != nullptr);
    args.insertOrAssign("material_names_path", DataPath({"Short Names"}));
    expectedCode = -12215;
  }
  const auto result = ComputeODFFilter().preflight(inputs.data_structure, args);
  REQUIRE(result.output_actions.invalid());
  REQUIRE(FirstErrorCode(result) == expectedCode);
}

TEST_CASE("MTRSim::ComputeODFFilter: Preflight defers crystal values for placeholder stores", "[MTRSim][ComputeODFFilter]")
{
  auto inputs = BuildSyntheticEbsd({{{11.0F, 22.0F, 33.0F}}}, {1}, k_HexagonalHigh);
  REQUIRE(inputs.data_structure.removeData(inputs.crystal_structures_path));
  const CreateArrayAction action(DataType::uint32, {2}, {1}, inputs.crystal_structures_path);
  const auto applied = action.apply(inputs.data_structure, IDataAction::Mode::Preflight);
  SIMPLNX_RESULT_REQUIRE_VALID(applied);
  REQUIRE(inputs.data_structure.removeData(inputs.material_names_path));
  const auto* ensemblePtr = inputs.data_structure.getDataAs<AttributeMatrix>(inputs.ensemble_attr_mat_path);
  REQUIRE(ensemblePtr != nullptr);
  auto* materialNamesPtr = StringArray::CreateWithValues(inputs.data_structure, "Material Names", {2}, {"", ""}, ensemblePtr->getId());
  REQUIRE(materialNamesPtr != nullptr);
  materialNamesPtr->setStore(std::make_shared<EmptyStringStore>(ShapeType{2}));
  const auto preflight = ComputeODFFilter().preflight(inputs.data_structure, MakeBaseArgs(inputs, false, 5.0F));
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.output_actions);
}

TEST_CASE("MTRSim::ComputeODFFilter: Create mode rejects the reserved metadata group name for cell data", "[MTRSim][ComputeODFFilter]")
{
  auto inputs = BuildSyntheticEbsd({{{11.0F, 22.0F, 33.0F}}}, {1}, k_HexagonalHigh);
  auto args = MakeBaseArgs(inputs, false, 5.0F);
  args.insertOrAssign(ComputeODFFilter::k_CellAttrMatName_Key, std::string{"ODF Metadata"});
  const auto preflight = ComputeODFFilter().preflight(inputs.data_structure, args);
  REQUIRE(preflight.output_actions.invalid());
  REQUIRE(FirstErrorCode(preflight) == -12218);
}

TEST_CASE("MTRSim::ComputeODFFilter: Execution validates selected crystal metadata before accumulation", "[MTRSim][ComputeODFFilter]")
{
  auto inputs = BuildSyntheticEbsd({{{11.0F, 22.0F, 33.0F}}}, {1}, 999);
  ComputeODFInputValues values{};
  values.eulerAnglesPath = inputs.euler_angles_path;
  values.phasesPath = inputs.phases_path;
  values.crystalStructuresPath = inputs.crystal_structures_path;
  int32 expectedCode = -12217;
  SECTION("unknown selected code")
  {
  }
  SECTION("placeholder crystal values")
  {
    REQUIRE(inputs.data_structure.removeData(inputs.crystal_structures_path));
    const CreateArrayAction action(DataType::uint32, {2}, {1}, inputs.crystal_structures_path);
    const auto applied = action.apply(inputs.data_structure, IDataAction::Mode::Preflight);
    SIMPLNX_RESULT_REQUIRE_VALID(applied);
  }
  SECTION("phase zero")
  {
    values.selectedPhase = 0;
    expectedCode = -12216;
  }
  SECTION("phase outside tuple range")
  {
    values.selectedPhase = 2;
    expectedCode = -12216;
  }
  const std::atomic_bool shouldCancel = false;
  const IFilter::MessageHandler handler;
  Result<> result;
  REQUIRE_NOTHROW(result = ComputeODF(inputs.data_structure, handler, shouldCancel, &values)());
  REQUIRE(result.invalid());
  REQUIRE(result.errors().front().code == expectedCode);
  REQUIRE(result.errors().front().message.find(inputs.crystal_structures_path.toString()) != std::string::npos);
}

TEST_CASE("MTRSim::ComputeODFFilter: Rejects metadata collisions", "[MTRSim][ComputeODFFilter]")
{
  auto inputs = BuildSyntheticEbsd({{{11.0F, 22.0F, 33.0F}}}, {1}, k_HexagonalHigh);
  BuildExistingOdfGeom(inputs.data_structure, {72, 36, 72}, {5.0F, 5.0F, 5.0F});
  const auto* geometryPtr = inputs.data_structure.getDataAs<ImageGeom>(k_OutputGeomPath);
  REQUIRE(geometryPtr != nullptr);
  int32 expectedCode = -12218;
  SECTION("wrong group type")
  {
    REQUIRE(StringArray::CreateWithValues(inputs.data_structure, "ODF Metadata", {1}, {"collision"}, geometryPtr->getId()) != nullptr);
  }
  SECTION("existing units object")
  {
    const auto* groupPtr = DataGroup::Create(inputs.data_structure, "ODF Metadata", geometryPtr->getId());
    REQUIRE(groupPtr != nullptr);
    REQUIRE(StringArray::CreateWithValues(inputs.data_structure, "Component 1 Units", {1}, {"MUD"}, groupPtr->getId()) != nullptr);
    expectedCode = -12219;
  }
  auto args = MakeBaseArgs(inputs, false, 5.0F);
  args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, ChoicesParameter::ValueType{1});
  args.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, k_OutputGeomPath);
  const auto result = ComputeODFFilter().preflight(inputs.data_structure, args);
  REQUIRE(result.output_actions.invalid());
  REQUIRE(FirstErrorCode(result) == expectedCode);
}

TEST_CASE("MTRSim::ComputeODFFilter: Units metadata survives DREAM3D round trip", "[MTRSim][ComputeODFFilter]")
{
  unit_test::LoadPlugins();
  const ChoicesParameter::ValueType units = GENERATE(0ULL, 1ULL);
  auto inputs = BuildSyntheticEbsd({{{11.0F, 22.0F, 33.0F}}}, {1}, k_HexagonalHigh);
  auto args = MakeBaseArgs(inputs, false, 10.0F);
  args.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, units);
  const auto computed = ComputeODFFilter().execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(computed.result);
  const std::filesystem::path outputPath = std::filesystem::path(unit_test::k_BinaryTestOutputDir.view()) / fmt::format("compute_odf_units_{}.dream3d", units);
  WriteDREAM3DFilter writer;
  auto writeArgs = writer.getDefaultArguments();
  writeArgs.insertOrAssign(WriteDREAM3DFilter::k_ExportFilePath, outputPath);
  writeArgs.insertOrAssign(WriteDREAM3DFilter::k_WriteXdmf, false);
  const auto written = writer.execute(inputs.data_structure, writeArgs);
  SIMPLNX_RESULT_REQUIRE_VALID(written.result);
  DataStructure restored;
  ReadDREAM3DFilter reader;
  auto readArgs = reader.getDefaultArguments();
  readArgs.insertOrAssign(ReadDREAM3DFilter::k_ImportFileData, Dream3dImportParameter::ImportData(outputPath));
  const auto read = reader.execute(restored, readArgs);
  SIMPLNX_RESULT_REQUIRE_VALID(read.result);
  const auto* unitsArrayPtr = restored.getDataAs<StringArray>(DataPath({"ODF", "ODF Metadata", "Component 1 Units"}));
  REQUIRE(unitsArrayPtr != nullptr);
  REQUIRE(unitsArrayPtr->getNumberOfTuples() == 1);
  REQUIRE((*unitsArrayPtr)[0] == (units == 0 ? "Count-Density" : "MUD"));
  REQUIRE(std::filesystem::remove(outputPath));
}

TEST_CASE("MTRSim::ComputeODFFilter: Append mode adds a new component to an "
          "existing ODF",
          "[MTRSim][ComputeODFFilter]")
{
  // Run Create New first: 1 HCP voxel at (10, 20, 30) deg, no smoothing, bin
  // size 5 deg. Expected Component 1 sum = 1.0 (12 deposits normalized by total
  // deposit count N=12, matches MATLAB calc_ODF.m).
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);
  {
    Arguments createArgs = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
    ComputeODFFilter filter;
    auto preflightResult = filter.preflight(inputs.data_structure, createArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
    auto executeResult = filter.execute(inputs.data_structure, createArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);
  }

  const DataPath component1Path = k_OutputGeomPath.createChildPath(k_CellAttrMatName).createChildPath(k_ComponentName);
  const auto [sum1Before, nonZero1Before] = SumAndNonZeroOf(inputs.data_structure, component1Path);
  REQUIRE(sum1Before == Approx(1.0).margin(1.0e-9));

  // Now Append a second component with the SAME EBSD input: sum should also
  // be 1.0.
  const std::string kAppendedComponentName = "Component 2";
  Arguments appendArgs = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  appendArgs.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  appendArgs.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, std::make_any<DataPath>(k_OutputGeomPath));
  appendArgs.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::make_any<DataObjectNameParameter::ValueType>(kAppendedComponentName));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, appendArgs);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, appendArgs);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const DataPath component2Path = k_OutputGeomPath.createChildPath(k_CellAttrMatName).createChildPath(kAppendedComponentName);
  REQUIRE_NOTHROW(inputs.data_structure.getDataRefAs<Float64ArrayType>(component2Path));
  const auto [sum2, nonZero2] = SumAndNonZeroOf(inputs.data_structure, component2Path);
  REQUIRE(sum2 == Approx(1.0).margin(1.0e-9));
  REQUIRE(nonZero2 == nonZero1Before);

  // Defense-in-depth: the appended array must have the same tuple count as
  // Component 1 (i.e. it was created against the existing geometry's dims, not
  // a fresh set).
  const auto& component1Arr = inputs.data_structure.getDataRefAs<Float64ArrayType>(component1Path);
  const auto& component2Arr = inputs.data_structure.getDataRefAs<Float64ArrayType>(component2Path);
  REQUIRE(component2Arr.getNumberOfTuples() == component1Arr.getNumberOfTuples());

  // Component 1 must be unchanged.
  const auto [sum1After, nonZero1After] = SumAndNonZeroOf(inputs.data_structure, component1Path);
  REQUIRE(sum1After == Approx(sum1Before).margin(1.0e-9));
  REQUIRE(nonZero1After == nonZero1Before);
}

TEST_CASE("MTRSim::ComputeODFFilter: Append mode rejects non-uniform spacing", "[MTRSim][ComputeODFFilter]")
{
  // Build EBSD inputs, then add a pre-existing ODF ImageGeom with non-uniform
  // spacing.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);
  BuildExistingOdfGeom(inputs.data_structure, {72, 36, 72}, {5.0f, 5.0f, 6.0f});

  Arguments appendArgs = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  appendArgs.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  appendArgs.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, std::make_any<DataPath>(k_OutputGeomPath));
  appendArgs.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::make_any<DataObjectNameParameter::ValueType>(std::string{"Component 2"}));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, appendArgs);
  REQUIRE(preflightResult.output_actions.invalid());
  bool found12207 = false;
  for(const auto& err : preflightResult.output_actions.errors())
  {
    if(err.code == -12207)
    {
      found12207 = true;
      // Also verify the error message is informative about the spacing issue.
      const bool messageMentionsSpacing =
          (err.message.find("uniform") != std::string::npos || err.message.find("non-uniform") != std::string::npos || err.message.find("spacing") != std::string::npos);
      REQUIRE(messageMentionsSpacing);
      break;
    }
  }
  REQUIRE(found12207);
}

TEST_CASE("MTRSim::ComputeODFFilter: Append mode rejects component-name collision", "[MTRSim][ComputeODFFilter]")
{
  // Existing ODF geom with a Component 1 seeded on it; Append call using the
  // same component name.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);
  BuildExistingOdfGeom(inputs.data_structure, {72, 36, 72}, {5.0f, 5.0f, 5.0f},
                       /*seedComponent1=*/true, /*seedValue=*/0.0);

  Arguments appendArgs = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  appendArgs.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  appendArgs.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, std::make_any<DataPath>(k_OutputGeomPath));
  // k_ComponentName_Key keeps its makeBaseArgs default = k_ComponentName =
  // "Component 1" → collision.

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, appendArgs);
  REQUIRE(preflightResult.output_actions.invalid());
  bool found12209 = false;
  for(const auto& err : preflightResult.output_actions.errors())
  {
    if(err.code == -12209)
    {
      found12209 = true;
      break;
    }
  }
  REQUIRE(found12209);
}

TEST_CASE("MTRSim::ComputeODFFilter: Append mode preflight reports derived bin size", "[MTRSim][ComputeODFFilter]")
{
  // Uniform 2.5-deg spacing existing geom; user passes bin_size_deg = 999.0 to
  // prove it's ignored.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);
  BuildExistingOdfGeom(inputs.data_structure, {144, 72, 144}, {2.5f, 2.5f, 2.5f});

  Arguments appendArgs = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/999.0f);
  appendArgs.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  appendArgs.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, std::make_any<DataPath>(k_OutputGeomPath));
  appendArgs.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::make_any<DataObjectNameParameter::ValueType>(std::string{"Component 2"}));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, appendArgs);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);

  const auto it = std::find_if(preflightResult.output_values.begin(), preflightResult.output_values.end(), [](const IFilter::PreflightValue& v) { return v.name == "Bin Size [deg]"; });
  REQUIRE(it != preflightResult.output_values.end());
  // Derived bin size should be 2.5, NOT 999.0.
  REQUIRE(it->value.find("2.5") != std::string::npos);
  REQUIRE(it->value.find("999") == std::string::npos);
}

TEST_CASE("MTRSim::ComputeODFFilter: Unknown selected crystal code is rejected", "[MTRSim][ComputeODFFilter]")
{
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, /*crystalCode=*/999u);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  REQUIRE(preflightResult.output_actions.invalid());
  REQUIRE(FirstErrorCode(preflightResult) == -12217);
  auto executeResult = filter.execute(inputs.data_structure, args);
  REQUIRE(executeResult.result.invalid());
}

TEST_CASE("MTRSim::ComputeODFFilter: Append mode rejects non-ImageGeom geometry", "[MTRSim][ComputeODFFilter]")
{
  // Build a DataStructure where /ODF resolves to a plain DataGroup (not an
  // ImageGeom), with a valid EBSD input tree alongside. Run Append mode
  // pointing at /ODF. The GeometrySelectionParameter wrapper enforces
  // IGeometry::Type::Image and flags this with its own generic geometry-type
  // code (-3) BEFORE preflightImpl runs, so our belt- and-suspenders -12206 in
  // preflightImpl is shadowed. Assert the wrapper code explicitly, mirroring
  // the pattern used by the "bad euler component count" test above.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);

  // Plant a DataGroup at /ODF instead of an ImageGeom.
  auto* dgPtr = DataGroup::Create(inputs.data_structure, "ODF");
  REQUIRE(dgPtr != nullptr);

  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key,
                      std::make_any<ChoicesParameter::ValueType>(1ULL)); // Append
  args.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  REQUIRE(preflightResult.output_actions.invalid());
  REQUIRE(FirstErrorCode(preflightResult) == -3);
}

TEST_CASE("MTRSim::ComputeODFFilter: Append mode rejects ImageGeom with no cell-data", "[MTRSim][ComputeODFFilter]")
{
  // Build a DataStructure with a bare ImageGeom at /ODF that has NO cell
  // AttributeMatrix attached. Run Append mode. Preflight must reject with code
  // -12208.
  EbsdInputs inputs = BuildSyntheticEbsd({{{10.0f, 20.0f, 30.0f}}}, {1}, k_HexagonalHigh);

  // Create a bare ImageGeom WITHOUT calling setCellData.
  auto* geomPtr = ImageGeom::Create(inputs.data_structure, "ODF");
  REQUIRE(geomPtr != nullptr);
  geomPtr->setDimensions({72, 36, 72});
  REQUIRE(geomPtr->setSpacing({5.0f, 5.0f, 5.0f}).valid());
  REQUIRE(geomPtr->setOrigin({0.0f, 0.0f, 0.0f}).valid());
  // Intentionally skip setCellData — leaves getCellData() == nullptr.

  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/false, /*binSizeDeg=*/5.0f);
  args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  args.insertOrAssign(ComputeODFFilter::k_ExistingOdfGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  REQUIRE(preflightResult.output_actions.invalid());
  REQUIRE(FirstErrorCode(preflightResult) == -12208);
}

// -----------------------------------------------------------------------------
// Phase 1: Targeted MATLAB calc_ODF.m bin-by-bin validation. The reference
// HDF5 is produced by matlab/run_validation.m which runs calc_ODF.m on the
// SAME 12-orientation table defined below; the test then loads it and diffs
// bin-by-bin against ComputeODFFilter's output. The reference HDF5 is NOT
// committed to the repo — the user runs MATLAB once locally to generate it.
// If the reference is missing, the test gracefully WARN-skips so CI passes.
// -----------------------------------------------------------------------------

namespace fs = std::filesystem;

TEST_CASE("MTRSim::ComputeODFFilter: Targeted MATLAB calc_ODF.m bin-by-bin "
          "validation",
          "[MTRSim][ComputeODFFilter]")
{
  // 12 hardcoded HCP orientations exercising PHI-near-pole, PHI-near-equator,
  // axis-wrap, and multi-axis generic cases. The reference HDF5 is produced
  // by matlab/run_validation.m which runs calc_ODF.m on the SAME 12-orientation
  // table defined below; the test then loads it and diffs bin-by-bin against
  // ComputeODFFilter's output.
  //
  // ---------------------------------------------------------------------------
  // FP-precision-safe fixture design (important):
  //   Every input value is of the form (5 deg * N + 2.5 deg) -- i.e. it sits
  //   at the EXACT MID-POINT of its 5 deg bin. Hex 6/mmm symmetry operators
  //   are at multiples of 60 deg (z-rot) and 30 deg (secondary 2-fold), so
  //   every symmetric variant of a mid-bin input is also mid-bin in all
  //   three Bunge axes. That gives 0.5 step-units (= 2.5 deg) of margin from
  //   the nearest bin boundary on every axis -- many orders of magnitude
  //   beyond the ulp-level FP precision of atan2/floor/fix at exact-boundary
  //   inputs. As a result, MATLAB / Python / C++ all agree bin-for-bin.
  //
  //   An earlier version of this fixture used inputs at exact bin boundaries
  //   (phi1=45, PHI=0/90/180, phi1/phi2=359, etc.). Those inputs produce
  //   variants whose Eulers fall within 1 ulp of bin boundaries, where IEEE
  //   754 atan2/cos/sin lib differences cause C++ and MATLAB to round the
  //   same mathematically-equivalent value to different bins -- ~1851 bins
  //   disagreed at max |diff| ~4e-3. Not an algorithm bug; pure FP precision
  //   pathology at the chosen inputs. Avoiding exact boundaries removes it.
  // ---------------------------------------------------------------------------
  //
  // This table is byte-identical to the comment + array in
  // matlab/run_validation.m. A future reader should be able to diff the two
  // and see they're the same.
  //
  //  # | (phi1 deg, PHI deg, phi2 deg) | bin (5 deg spacing)   | Why
  // ---+--------------------------------+-----------------------+----------------------------------------------
  //   1| ( 12.5,    12.5,   12.5)       | (2,  2,  2)           | Pure mid-bin;
  //   baseline 2| ( 47.5,    27.5,   92.5)       | (9,  5,  18)          |
  //   Generic interior, no special structure 3| (137.5,    67.5,  217.5) | (27,
  //   13, 43)          | Multi-decimal-bin coverage; phi1/phi2 > 90 deg 4|
  //   (  2.5,     2.5,    2.5)       | (0,  0,  0)           | All near zero;
  //   tests near-PHI=0 handling 5| ( 47.5,     2.5,   32.5)       | (9,  0,  6)
  //   | PHI very small with non-trivial phi1/phi2 6| ( 62.5,     2.5,   92.5)
  //   | (12, 0,  18)          | PHI very small, generic phi1/phi2 7|
  //   ( 47.5,    92.5,   32.5)       | (9,  18, 6)           | PHI just above
  //   pi/2 (equatorial regime) 8| ( 47.5,   177.5,   32.5)       | (9,  35, 6)
  //   | PHI near pi (upper-pole regime) 9| ( 47.5,    87.5,   32.5)       | (9,
  //   17, 6)           | PHI just below pi/2 (complement of #7)
  //  10| (357.5,    32.5,   32.5)       | (71, 6,  6)           | phi1 near 360
  //  deg (wrap regime) 11| ( 47.5,    32.5,  357.5)       | (9,  6,  71) | phi2
  //  near 360 deg (wrap regime) 12| (357.5,    87.5,  357.5)       | (71, 17,
  //  71)          | All three near upper boundaries simultaneously
  const std::vector<std::array<float32, 3>> eulersDeg = {
      {12.5f, 12.5f, 12.5f}, {47.5f, 27.5f, 92.5f},  {137.5f, 67.5f, 217.5f}, {2.5f, 2.5f, 2.5f},     {47.5f, 2.5f, 32.5f},   {62.5f, 2.5f, 92.5f},
      {47.5f, 92.5f, 32.5f}, {47.5f, 177.5f, 32.5f}, {47.5f, 87.5f, 32.5f},   {357.5f, 32.5f, 32.5f}, {47.5f, 32.5f, 357.5f}, {357.5f, 87.5f, 357.5f},
  };

  // buildSyntheticEbsd already converts deg -> rad internally (see helper).
  // It mirrors what calc_ODF.m receives (phi1/PHI/phi2 in radians).
  EbsdInputs inputs = BuildSyntheticEbsd(eulersDeg, std::vector<int32>(eulersDeg.size(), 1), k_HexagonalHigh);
  Arguments args = MakeBaseArgs(inputs, /*applySmoothing=*/true, /*binSizeDeg=*/5.0f);

  ComputeODFFilter filter;
  auto preflightResult = filter.preflight(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(inputs.data_structure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // Load the MATLAB reference. The reference HDF5 IS committed (see commit
  // 03ed497) but is regenerated whenever the orientation table above changes.
  //
  // Two skip-with-WARN paths:
  //   1. File missing entirely -> user never ran MATLAB locally.
  //   2. File present but its embedded fixture_version doesn't match
  //      k_TargetedFixtureVersion -> user ran MATLAB on an older fixture and
  //      the reference is now stale.
  // Both cases prompt the user to (re)run run_validation.m and regenerate.
  //
  // Bumping rules: when changing the orientation table, bump BOTH:
  //   - This constant
  //   - The matching `fixture_version` write in matlab/run_validation.m
  //   v1: original boundary-stress fixture (deprecated; FP-precision pathology)
  //   v2: FP-precision-safe mid-bin fixture
  constexpr int64 k_TargetedFixtureVersion = 2;
  const fs::path refPath = fs::path(fmt::format("{}/calc_odf_reference_targeted.h5", unit_test::k_DataDir.view()));
  const std::string regenerateMsg = "run `matlab -batch \"addpath('matlab'); run_validation();\"` "
                                    "from the repo root to (re)generate it. Skipping bin-by-bin comparison.";
  if(!fs::exists(refPath))
  {
    WARN("Reference HDF5 not present at " << refPath.string() << " - " << regenerateMsg);
    return;
  }
  const std::optional<int64> refFixtureVersion = mtrsim::tryReadFixtureVersion(refPath.string());
  if(!refFixtureVersion.has_value())
  {
    WARN("Reference HDF5 at " << refPath.string() << " has no /ODF_best/fixture_version field "
                              << "(produced by an older run_validation.m "
                                 "before the fixture-version mechanism). "
                              << regenerateMsg);
    return;
  }
  if(*refFixtureVersion != k_TargetedFixtureVersion)
  {
    WARN("Reference HDF5 fixture_version mismatch: file has v" << *refFixtureVersion << " but test expects v" << k_TargetedFixtureVersion << ". " << regenerateMsg);
    return;
  }

  const auto refComponents = mtrsim::readODFComponents(refPath.string());
  REQUIRE(refComponents.size() == 1);
  const auto& refValues = refComponents[0].values;

  // Compute output ODF array (default geometry/cell-attr-mat/component path,
  // mirroring sumAndNonZeroCount above).
  const DataPath odfPath = k_OutputGeomPath.createChildPath(k_CellAttrMatName).createChildPath(k_ComponentName);
  const auto& outArr = inputs.data_structure.getDataRefAs<Float64ArrayType>(odfPath);
  const auto& outStore = outArr.getDataStoreRef();
  REQUIRE(outStore.getSize() == refValues.size());

  // Bin-by-bin diff. Tolerance 1e-9 -- both sides do the same arithmetic in
  // double precision, so any difference > 1e-9 means a real bin mismatch
  // (PHI clamp vs wrap, an off-by-one in the smoothing kernel, etc.) rather
  // than floating-point noise.
  double maxAbsDiff = 0.0;
  double sumSqDiff = 0.0;
  std::size_t failingBins = 0;
  std::size_t firstFailingIndex = SIZE_MAX;
  constexpr double k_Tolerance = 1.0e-9;
  for(std::size_t i = 0; i < refValues.size(); ++i)
  {
    const double diff = std::abs(static_cast<double>(outStore[i]) - refValues[i]);
    maxAbsDiff = std::max(diff, maxAbsDiff);
    sumSqDiff += diff * diff;
    if(diff > k_Tolerance)
    {
      if(firstFailingIndex == SIZE_MAX)
      {
        firstFailingIndex = i;
      }
      ++failingBins;
    }
  }
  const double rmsDiff = std::sqrt(sumSqDiff / static_cast<double>(refValues.size()));

  INFO(fmt::format("max |diff| = {:.3e}, RMS = {:.3e}, failing bins = {} / {} "
                   "(tolerance = {:.0e})",
                   maxAbsDiff, rmsDiff, failingBins, refValues.size(), k_Tolerance));
  if(failingBins > 0)
  {
    INFO(fmt::format("First failing bin: index {}, MATLAB = {:.6e}, C++ = "
                     "{:.6e}, diff = {:.3e}",
                     firstFailingIndex, refValues[firstFailingIndex], static_cast<double>(outStore[firstFailingIndex]),
                     std::abs(static_cast<double>(outStore[firstFailingIndex]) - refValues[firstFailingIndex])));
  }
  REQUIRE(failingBins == 0);
}

TEST_CASE("MTRSim::ComputeODFFilter: Realistic 640x640 HCP scan validation "
          "against calc_ODF.m",
          "[MTRSim][ComputeODFFilter][validation]")
{
  // Phase 2 validation: load the 640x640 HCP titanium EBSD scan from
  // data/real_world_microtexture_data.dream3d, run ComputeODFFilter on the
  // ~357k masked phase-1 voxels, and bin-by-bin diff against the MATLAB
  // calc_ODF.m reference at data/calc_odf_reference_realistic.h5.
  //
  // Tagged [validation] so it can be skipped from the fast inner-loop test
  // suite via `ctest -E validation`. Run this explicitly with
  // `ctest -R "Realistic" --output-on-failure`.
  //
  // The .dream3d file's Mask is stored as uint8 (DREAM3D-NX convention) but
  // ComputeODFFilter expects a Bool mask, so we build a Bool copy in-place
  // before invoking the filter.
  unit_test::LoadPlugins();

  const fs::path inputFile = fs::path(fmt::format("{}/real_world_microtexture_data.dream3d", unit_test::k_DataDir.view()));
  const fs::path refFile = fs::path(fmt::format("{}/calc_odf_reference_realistic.h5", unit_test::k_DataDir.view()));
  if(!fs::exists(inputFile))
  {
    WARN("Realistic EBSD fixture not present at " << inputFile.string() << " - skipping.");
    return;
  }
  if(!fs::exists(refFile))
  {
    WARN("Realistic MATLAB reference not present at " << refFile.string() << " - run matlab/run_validation.m to generate it. Skipping.");
    return;
  }

  // --- Load the .dream3d via ReadDREAM3DFilter
  DataStructure dataStructure;
  {
    ReadDREAM3DFilter readFilter;
    Arguments readArgs;
    Dream3dImportParameter::ImportData importData(inputFile);
    readArgs.insertOrAssign(ReadDREAM3DFilter::k_ImportFileData, std::make_any<Dream3dImportParameter::ImportData>(importData));
    auto preflight = readFilter.preflight(dataStructure, readArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(preflight.output_actions);
    auto execute = readFilter.execute(dataStructure, readArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(execute.result);
  }

  // Paths inside the loaded DataStructure (per the .dream3d's internal layout).
  const DataPath geomPath({"16_Micro_Res"});
  const DataPath cellAttrMatPath = geomPath.createChildPath("CellData");
  const DataPath eulerAnglesPath = cellAttrMatPath.createChildPath("EulerAngles");
  const DataPath phasesPath = cellAttrMatPath.createChildPath("Phases");
  const DataPath maskUInt8Path = cellAttrMatPath.createChildPath("Mask");
  const DataPath crystalStructuresPath = geomPath.createChildPath("CellEnsembleData").createChildPath("CrystalStructures");
  const DataPath maskBoolPath = cellAttrMatPath.createChildPath("MaskBool");

  // Build a Bool mask alongside whatever-type Mask was loaded.
  // ReadDREAM3DFilter may give it back as UInt8Array or BoolArray depending on
  // stored metadata — handle both so we don't bad_cast on either.
  {
    const auto* maskBasePtr = dataStructure.getDataAs<IDataArray>(maskUInt8Path);
    REQUIRE(maskBasePtr != nullptr);
    INFO("Loaded Mask type name: " << maskBasePtr->getTypeName());

    const usize numTuples = maskBasePtr->getNumberOfTuples();
    auto* cellAmPtr = dataStructure.getDataAs<AttributeMatrix>(cellAttrMatPath);
    REQUIRE(cellAmPtr != nullptr);
    auto* maskBoolPtr = BoolArrayType::CreateWithStore<DataStore<bool>>(dataStructure, "MaskBool", {numTuples}, {1}, cellAmPtr->getId());
    REQUIRE(maskBoolPtr != nullptr);
    auto& boolStore = maskBoolPtr->getDataStoreRef();

    if(const auto* maskU8Ptr = dataStructure.getDataAs<UInt8ArrayType>(maskUInt8Path); maskU8Ptr != nullptr)
    {
      const auto& u8Store = maskU8Ptr->getDataStoreRef();
      for(usize i = 0; i < numTuples; ++i)
      {
        boolStore.setValue(i, u8Store.getValue(i) != 0);
      }
    }
    else if(const auto* maskBlPtr = dataStructure.getDataAs<BoolArrayType>(maskUInt8Path); maskBlPtr != nullptr)
    {
      const auto& blStore = maskBlPtr->getDataStoreRef();
      for(usize i = 0; i < numTuples; ++i)
      {
        boolStore.setValue(i, blStore.getValue(i));
      }
    }
    else
    {
      FAIL("Mask array at " + maskUInt8Path.toString() + " is neither UInt8Array nor BoolArray");
    }
  }

  // --- Run ComputeODFFilter
  ComputeODFFilter filter;
  Arguments args;
  args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  args.insertOrAssign(ComputeODFFilter::k_ApplySmoothing_Key, std::make_any<bool>(true));
  args.insertOrAssign(ComputeODFFilter::k_BinSizeDeg_Key, std::make_any<float32>(5.0f));
  args.insertOrAssign(ComputeODFFilter::k_EulerAngles_Key, std::make_any<DataPath>(eulerAnglesPath));
  args.insertOrAssign(ComputeODFFilter::k_Phases_Key, std::make_any<DataPath>(phasesPath));
  args.insertOrAssign(ComputeODFFilter::k_CrystalStructures_Key, std::make_any<DataPath>(crystalStructuresPath));
  args.insertOrAssign(ComputeODFFilter::k_MaterialNames_Key, geomPath.createChildPath("CellEnsembleData").createChildPath("MaterialName"));
  args.insertOrAssign(ComputeODFFilter::k_PhaseIndex_Key, int32{1});
  args.insertOrAssign(ComputeODFFilter::k_UseMask_Key, std::make_any<bool>(true));
  args.insertOrAssign(ComputeODFFilter::k_MaskArrayPath_Key, std::make_any<DataPath>(maskBoolPath));
  args.insertOrAssign(ComputeODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(ComputeODFFilter::k_CellAttrMatName_Key, std::make_any<DataObjectNameParameter::ValueType>(std::string{"Cell Data"}));
  args.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::make_any<DataObjectNameParameter::ValueType>(std::string{"Component 1"}));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // --- Diff bin-by-bin against the MATLAB reference
  const auto refComponents = mtrsim::readODFComponents(refFile);
  REQUIRE(refComponents.size() == 1);
  const auto& refValues = refComponents[0].values;

  const auto& outArr = dataStructure.getDataRefAs<Float64ArrayType>(DataPath({"ODF", "Cell Data", "Component 1"}));
  const auto& outStore = outArr.getDataStoreRef();
  REQUIRE(outStore.getSize() == refValues.size());

  double maxAbsDiff = 0.0;
  double sumSqDiff = 0.0;
  std::size_t failingBins = 0;
  std::size_t firstFailingIndex = SIZE_MAX;
  constexpr double k_Tolerance = 1.0e-9;
  for(std::size_t i = 0; i < refValues.size(); ++i)
  {
    const double diff = std::abs(static_cast<double>(outStore[i]) - refValues[i]);
    maxAbsDiff = std::max(diff, maxAbsDiff);
    sumSqDiff += diff * diff;
    if(diff > k_Tolerance)
    {
      if(firstFailingIndex == SIZE_MAX)
      {
        firstFailingIndex = i;
      }
      ++failingBins;
    }
  }
  const double rmsDiff = std::sqrt(sumSqDiff / static_cast<double>(refValues.size()));

  INFO(fmt::format("max |diff| = {:.3e}, RMS = {:.3e}, failing bins = {} / {} "
                   "(tolerance = {:.0e})",
                   maxAbsDiff, rmsDiff, failingBins, refValues.size(), k_Tolerance));
  if(failingBins > 0)
  {
    INFO(fmt::format("First failing bin: index {}, MATLAB = {:.6e}, C++ = "
                     "{:.6e}, diff = {:.3e}",
                     firstFailingIndex, refValues[firstFailingIndex], static_cast<double>(outStore[firstFailingIndex]),
                     std::abs(static_cast<double>(outStore[firstFailingIndex]) - refValues[firstFailingIndex])));
  }
  REQUIRE(failingBins == 0);
}

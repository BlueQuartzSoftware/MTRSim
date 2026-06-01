/**
 * Unit tests for MTRSimFilter (Milestone AJ, Tasks 5 + 6 scaffold).
 *
 * Covers preflight parameter validation only (the algorithm body is a no-op
 * stub at this stage):
 *   1. Volume Fraction column count != numComponents -> invalid.
 *   2. Theta List rows < numComponents - 1 -> invalid.
 *   3. Volume Fraction values do not sum to 1.0 -> invalid.
 *   4. Happy-path args -> preflight VALID (builds geometry + array actions).
 *   5. Theta List rows with wrong column count (2 instead of 3) -> invalid (-13005).
 */

#include <catch2/catch.hpp>

#include "MTRSim/Filters/MTRSimFilter.hpp"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/DynamicTableParameter.hpp"
#include "simplnx/Parameters/MultiArraySelectionParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <fmt/format.h>

#include <array>
#include <cmath>
#include <vector>

using namespace nx::core;
using namespace nx::core::UnitTest;

namespace
{
constexpr usize k_NPhi1 = 72;
constexpr usize k_NPHI = 36;
constexpr usize k_NPhi2 = 72;

const DataPath k_OdfGeomPath({"ODF"});
const std::string k_CellAttrMatName = "Cell Data";

// Builds an ODF ImageGeom (72x36x72, 5-degree spacing) with numComponents
// Float64 single-component cell arrays named component_0.., returning their
// DataPaths.
std::vector<DataPath> BuildOdfDataStructure(DataStructure& dataStructure, usize numComponents)
{
  ImageGeom* imageGeom = ImageGeom::Create(dataStructure, k_OdfGeomPath.getTargetName());
  imageGeom->setSpacing({5.0f, 5.0f, 5.0f});
  imageGeom->setOrigin({0.0f, 0.0f, 0.0f});
  imageGeom->setDimensions({k_NPhi2, k_NPHI, k_NPhi1}); // X(phi2), Y(PHI), Z(phi1)

  // ZYX tuple shape (slowest to fastest)
  const ShapeType tupleShape = {k_NPhi1, k_NPHI, k_NPhi2};

  AttributeMatrix* cellAM = AttributeMatrix::Create(dataStructure, k_CellAttrMatName, tupleShape, imageGeom->getId());

  std::vector<DataPath> compPaths;
  const DataPath cellAttrMatPath = k_OdfGeomPath.createChildPath(k_CellAttrMatName);
  for(usize c = 0; c < numComponents; ++c)
  {
    const std::string name = fmt::format("component_{}", c);
    CreateTestDataArray<float64>(dataStructure, name, tupleShape, {1}, cellAM->getId());
    compPaths.push_back(cellAttrMatPath.createChildPath(name));
  }
  return compPaths;
}

// Builds a valid argument set for the supplied component paths.
Arguments MakeValidArgs(const std::vector<DataPath>& compPaths)
{
  Arguments args;
  args.insertOrAssign(MTRSimFilter::k_InputOdfGeometry_Key, k_OdfGeomPath);
  args.insertOrAssign(MTRSimFilter::k_OdfComponentArrays_Key, compPaths);
  args.insertOrAssign(MTRSimFilter::k_VolumeFractions_Key, DynamicTableParameter::ValueType{{0.30, 0.35, 0.35}});
  args.insertOrAssign(MTRSimFilter::k_ThetaList_Key, DynamicTableParameter::ValueType{{0.1, 0.45, 0.1}, {0.08, 0.37, 0.08}});
  args.insertOrAssign(MTRSimFilter::k_PhysicalSize_Key, std::vector<float32>{2.0f, 2.0f, 0.0f});
  args.insertOrAssign(MTRSimFilter::k_PhysicalSpacing_Key, std::vector<float32>{0.02f, 0.02f, 0.02f});
  args.insertOrAssign(MTRSimFilter::k_UseSeed_Key, true);
  args.insertOrAssign(MTRSimFilter::k_SeedValue_Key, static_cast<uint64>(42));
  args.insertOrAssign(MTRSimFilter::k_SeedArrayName_Key, std::string("MTRSim SeedValue"));
  args.insertOrAssign(MTRSimFilter::k_GeneratePolarColoring_Key, false);
  args.insertOrAssign(MTRSimFilter::k_OutputGeometry_Key, DataPath({"MTR Microstructure"}));
  args.insertOrAssign(MTRSimFilter::k_CellAttrMatName_Key, std::string("Cell Data"));
  args.insertOrAssign(MTRSimFilter::k_MtrIdsArrayName_Key, std::string("MTRIds"));
  args.insertOrAssign(MTRSimFilter::k_EulersArrayName_Key, std::string("Eulers"));
  args.insertOrAssign(MTRSimFilter::k_PolarColorsArrayName_Key, std::string("Polar Colors"));
  return args;
}
} // namespace

TEST_CASE("MTRSim::MTRSimFilter: Valid preflight builds actions", "[MTRSim][MTRSimFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;
  const std::vector<DataPath> compPaths = BuildOdfDataStructure(dataStructure, 3);

  MTRSimFilter filter;
  Arguments args = MakeValidArgs(compPaths);

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
}

TEST_CASE("MTRSim::MTRSimFilter: Rejects mismatched Volume Fraction column count", "[MTRSim][MTRSimFilter][ErrorPath]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;
  const std::vector<DataPath> compPaths = BuildOdfDataStructure(dataStructure, 3);

  MTRSimFilter filter;
  Arguments args = MakeValidArgs(compPaths);
  // Only 2 volume-fraction columns for 3 components.
  args.insertOrAssign(MTRSimFilter::k_VolumeFractions_Key, DynamicTableParameter::ValueType{{0.5, 0.5}});

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.outputActions);
}

TEST_CASE("MTRSim::MTRSimFilter: Execute wires simulation to output arrays", "[MTRSim][MTRSimFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;
  const std::vector<DataPath> compPaths = BuildOdfDataStructure(dataStructure, 3);
  // Fill every component array with a uniform value so ODF sampling is well-defined.
  for(const auto& path : compPaths)
  {
    auto& arr = dataStructure.getDataRefAs<Float64Array>(path);
    arr.fill(1.0);
  }

  MTRSimFilter filter;
  Arguments args = MakeValidArgs(compPaths);
  // size/spacing/VF come from MakeValidArgs (100x100 grid)
  args.insertOrAssign(MTRSimFilter::k_UseSeed_Key, true);
  args.insertOrAssign(MTRSimFilter::k_SeedValue_Key, static_cast<uint64>(42));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // Output ImageGeom exists.
  const DataPath outGeomPath({"MTR Microstructure"});
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(outGeomPath));

  const DataPath cellAm = outGeomPath.createChildPath("Cell Data");
  const usize expectedTuples = 100 * 100;

  // MTRIds: Int32, 1 component, 10000 tuples.
  auto& mtrIds = dataStructure.getDataRefAs<Int32Array>(cellAm.createChildPath("MTRIds"));
  REQUIRE(mtrIds.getNumberOfComponents() == 1);
  REQUIRE(mtrIds.getNumberOfTuples() == expectedTuples);

  // Eulers: Float32, 3 components, 10000 tuples.
  auto& eulers = dataStructure.getDataRefAs<Float32Array>(cellAm.createChildPath("Eulers"));
  REQUIRE(eulers.getNumberOfComponents() == 3);
  REQUIRE(eulers.getNumberOfTuples() == expectedTuples);

  // MTR ids in {1,2,3}; at least 2 distinct ids appear. Also accumulate empirical
  // volume fractions for a loose wiring check.
  const auto& mtrStore = mtrIds.getDataStoreRef();
  std::array<usize, 4> counts = {0, 0, 0, 0};
  for(usize i = 0; i < mtrStore.getSize(); ++i)
  {
    const int32 id = mtrStore[i];
    REQUIRE(id >= 1);
    REQUIRE(id <= 3);
    counts[static_cast<usize>(id)]++;
  }
  usize distinct = 0;
  for(usize id = 1; id <= 3; ++id)
  {
    if(counts[id] > 0)
    {
      distinct++;
    }
  }
  REQUIRE(distinct >= 2);

  // Euler values finite and within Bunge bounds (interleaved 3/voxel).
  constexpr float twoPi = 2.0f * static_cast<float>(M_PI);
  constexpr float pi = static_cast<float>(M_PI);
  const auto& eulerStore = eulers.getDataStoreRef();
  for(usize t = 0; t < expectedTuples; ++t)
  {
    const float phi1 = eulerStore[t * 3 + 0];
    const float Phi = eulerStore[t * 3 + 1];
    const float phi2 = eulerStore[t * 3 + 2];
    REQUIRE(std::isfinite(phi1));
    REQUIRE(std::isfinite(Phi));
    REQUIRE(std::isfinite(phi2));
    REQUIRE(phi1 >= 0.0f);
    REQUIRE(phi1 <= twoPi);
    REQUIRE(Phi >= 0.0f);
    REQUIRE(Phi <= pi);
    REQUIRE(phi2 >= 0.0f);
    REQUIRE(phi2 <= twoPi);
  }

  // Seed array records 42.
  auto& seedArray = dataStructure.getDataRefAs<UInt64Array>(DataPath({"MTRSim SeedValue"}));
  REQUIRE(seedArray[0] == 42);

  // Loose volume-fraction wiring check (NOT a statistics check; rigorous VF
  // validation lives in the LibMTRSim statistical test). A 100x100 correlated
  // field has real variance, so use a generous margin of 0.12.
  const std::array<double, 4> targets = {0.0, 0.30, 0.35, 0.35};
  for(usize id = 1; id <= 3; ++id)
  {
    const double empirical = static_cast<double>(counts[id]) / static_cast<double>(expectedTuples);
    REQUIRE(empirical == Approx(targets[id]).margin(0.12));
  }
}

TEST_CASE("MTRSim::MTRSimFilter: Rejects too few Theta List rows", "[MTRSim][MTRSimFilter][ErrorPath]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;
  const std::vector<DataPath> compPaths = BuildOdfDataStructure(dataStructure, 3);

  MTRSimFilter filter;
  Arguments args = MakeValidArgs(compPaths);
  // 3 components require >= 2 theta rows; supply only 1.
  args.insertOrAssign(MTRSimFilter::k_ThetaList_Key, DynamicTableParameter::ValueType{{0.1, 0.45, 0.1}});

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.outputActions);
}

TEST_CASE("MTRSim::MTRSimFilter: Rejects Volume Fraction not summing to 1.0", "[MTRSim][MTRSimFilter][ErrorPath]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;
  const std::vector<DataPath> compPaths = BuildOdfDataStructure(dataStructure, 3);

  MTRSimFilter filter;
  Arguments args = MakeValidArgs(compPaths);
  // Columns sum to 0.6, not 1.0.
  args.insertOrAssign(MTRSimFilter::k_VolumeFractions_Key, DynamicTableParameter::ValueType{{0.2, 0.2, 0.2}});

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.outputActions);
}

TEST_CASE("MTRSim::MTRSimFilter: Rejects Theta List rows with wrong column count", "[MTRSim][MTRSimFilter][ErrorPath]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;
  const std::vector<DataPath> compPaths = BuildOdfDataStructure(dataStructure, 3);

  MTRSimFilter filter;
  Arguments args = MakeValidArgs(compPaths);
  // 2 rows supplied (enough for 3 components: needs >= 2), but each row has only 2
  // columns instead of 3 — must trigger the column-count check (-13005), not the
  // row-count check (-13004).
  args.insertOrAssign(MTRSimFilter::k_ThetaList_Key, DynamicTableParameter::ValueType{{0.1, 0.45}, {0.08, 0.37}});

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.outputActions);
}

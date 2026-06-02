/**
 * Unit tests for ReadMTRSimODFFilter (Milestone AJ, Task 3 + Task 2 path-prefix
 * cross-cut).
 *
 * Covers:
 *   1. Happy path on the canonical /ODF_best exemplar, with assertions on the
 *      new PreflightUpdatedValues ("HDF5 Path Prefix" label) and explicit
 * prefix arg.
 *   2. Error path: non-existent file -> preflight returns invalid.
 *   3. Blank ODF fixture (prefix "/blank_ODF"): execute succeeds, one
 * zero-valued Float64 component array is created.
 *   4. Uniform ODF fixture (prefix "/uniform_ODF"): execute succeeds, component
 *      sums to ~1.0 and every value is strictly positive (the fixture is a
 *      normalised-over-FZ ODF, NOT a per-cell uniform distribution — see repo
 *      TODO comments for rationale).
 *   5. Error path: wrong prefix against simulation_ODF.h5 -> preflight invalid.
 */

#include <catch2/catch.hpp>

#include "MTRSim/Filters/ReadMTRSimODFFilter.hpp"
#include "MTRSim/MTRSim_test_dirs.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <filesystem>

namespace fs = std::filesystem;

using namespace nx::core;
using namespace nx::core::UnitTest;

namespace
{
constexpr usize k_ExpectedNPhi1 = 72;
constexpr usize k_ExpectedNPHI = 36;
constexpr usize k_ExpectedNPhi2 = 72;
constexpr usize k_ExpectedNumTuples = k_ExpectedNPhi1 * k_ExpectedNPHI * k_ExpectedNPhi2; // 186624
constexpr float32 k_ExpectedSpacingDeg = 5.0f;
constexpr int32 k_ExpectedNumComponents = 3;
} // namespace

TEST_CASE("MTRSim::ReadMTRSimODFFilter: Valid Filter Execution", "[MTRSim][ReadMTRSimODFFilter]")
{
  UnitTest::LoadPlugins();

  const fs::path inputFile = fs::path(fmt::format("{}/simulation_ODF.h5", unit_test::k_DataDir.view()));
  REQUIRE(fs::exists(inputFile));

  const DataPath imageGeomPath({"ODF"});
  const std::string cellAttrMatName = "Cell Data";

  DataStructure dataStructure;

  ReadMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(inputFile));
  args.insertOrAssign(ReadMTRSimODFFilter::k_Hdf5PathPrefix_Key, std::make_any<StringParameter::ValueType>(std::string{"/ODF_best"}));
  args.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(imageGeomPath));
  args.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key, std::make_any<DataObjectNameParameter::ValueType>(cellAttrMatName));

  // Preflight
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  // PreflightUpdatedValues: must surface the prefix label so the UI can preview
  // it.
  REQUIRE(!preflightResult.outputValues.empty());
  const bool foundPrefixLabel = std::any_of(preflightResult.outputValues.begin(), preflightResult.outputValues.end(), [](const IFilter::PreflightValue& v) { return v.name == "HDF5 Path Prefix"; });
  REQUIRE(foundPrefixLabel);

  // Execute
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // Verify the ImageGeom was created at /ODF with the expected dims & spacing.
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(imageGeomPath));
  const auto& imageGeom = dataStructure.getDataRefAs<ImageGeom>(imageGeomPath);

  const auto dims = imageGeom.getDimensions();
  REQUIRE(dims[0] == k_ExpectedNPhi2); // X <- phi2
  REQUIRE(dims[1] == k_ExpectedNPHI);  // Y <- PHI
  REQUIRE(dims[2] == k_ExpectedNPhi1); // Z <- phi1

  const auto spacing = imageGeom.getSpacing();
  REQUIRE(std::fabs(spacing[0] - k_ExpectedSpacingDeg) < 1e-6f);
  REQUIRE(std::fabs(spacing[1] - k_ExpectedSpacingDeg) < 1e-6f);
  REQUIRE(std::fabs(spacing[2] - k_ExpectedSpacingDeg) < 1e-6f);

  // Verify each per-component Float64 array exists with the expected tuple
  // count and that component_0's values sum to ~1.0 (normalized ODF, within the
  // 1% MATLAB tolerance).
  const DataPath cellAttrMatPath = imageGeomPath.createChildPath(cellAttrMatName);
  for(int32 c = 0; c < k_ExpectedNumComponents; ++c)
  {
    const DataPath componentPath = cellAttrMatPath.createChildPath(fmt::format("component_{}", c));
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(componentPath));
    const auto& componentArray = dataStructure.getDataRefAs<Float64Array>(componentPath);
    REQUIRE(componentArray.getNumberOfTuples() == k_ExpectedNumTuples);
    REQUIRE(componentArray.getNumberOfComponents() == 1);
  }

  const auto& component0 = dataStructure.getDataRefAs<Float64Array>(cellAttrMatPath.createChildPath("component_0"));
  double sum0 = 0.0;
  for(usize i = 0; i < component0.getSize(); ++i)
  {
    sum0 += component0[i];
  }
  REQUIRE(std::fabs(sum0 - 1.0) < 1e-2);
}

TEST_CASE("MTRSim::ReadMTRSimODFFilter: Missing File", "[MTRSim][ReadMTRSimODFFilter][ErrorPath]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;

  ReadMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(fs::path("/nonexistent.h5")));
  args.insertOrAssign(ReadMTRSimODFFilter::k_Hdf5PathPrefix_Key, std::make_any<StringParameter::ValueType>(std::string{"/ODF_best"}));
  args.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key, std::make_any<DataObjectNameParameter::ValueType>(std::string{"Cell Data"}));

  auto preflightResult = filter.preflight(dataStructure, args);
  REQUIRE(preflightResult.outputActions.invalid());
}

TEST_CASE("MTRSim::ReadMTRSimODFFilter: Reads blank ODF fixture", "[MTRSim][ReadMTRSimODFFilter]")
{
  UnitTest::LoadPlugins();

  const fs::path inputFile = fs::path(fmt::format("{}/blank_ODF.h5", unit_test::k_DataDir.view()));
  REQUIRE(fs::exists(inputFile));

  const DataPath imageGeomPath({"ODF"});
  const std::string cellAttrMatName = "Cell Data";

  DataStructure dataStructure;

  ReadMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(inputFile));
  args.insertOrAssign(ReadMTRSimODFFilter::k_Hdf5PathPrefix_Key, std::make_any<StringParameter::ValueType>(std::string{"/blank_ODF"}));
  args.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(imageGeomPath));
  args.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key, std::make_any<DataObjectNameParameter::ValueType>(cellAttrMatName));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(imageGeomPath));
  const auto& imageGeom = dataStructure.getDataRefAs<ImageGeom>(imageGeomPath);
  const auto dims = imageGeom.getDimensions();
  REQUIRE(dims[0] == k_ExpectedNPhi2);
  REQUIRE(dims[1] == k_ExpectedNPHI);
  REQUIRE(dims[2] == k_ExpectedNPhi1);

  const DataPath componentPath = imageGeomPath.createChildPath(cellAttrMatName).createChildPath("component_0");
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(componentPath));
  const auto& component0 = dataStructure.getDataRefAs<Float64Array>(componentPath);
  REQUIRE(component0.getNumberOfTuples() == k_ExpectedNumTuples);

  for(usize i = 0; i < component0.getSize(); ++i)
  {
    REQUIRE(component0[i] == 0.0);
  }
}

TEST_CASE("MTRSim::ReadMTRSimODFFilter: Reads uniform ODF fixture", "[MTRSim][ReadMTRSimODFFilter]")
{
  UnitTest::LoadPlugins();

  const fs::path inputFile = fs::path(fmt::format("{}/uniform_ODF.h5", unit_test::k_DataDir.view()));
  REQUIRE(fs::exists(inputFile));

  const DataPath imageGeomPath({"ODF"});
  const std::string cellAttrMatName = "Cell Data";

  DataStructure dataStructure;

  ReadMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(inputFile));
  args.insertOrAssign(ReadMTRSimODFFilter::k_Hdf5PathPrefix_Key, std::make_any<StringParameter::ValueType>(std::string{"/uniform_ODF"}));
  args.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(imageGeomPath));
  args.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key, std::make_any<DataObjectNameParameter::ValueType>(cellAttrMatName));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const DataPath componentPath = imageGeomPath.createChildPath(cellAttrMatName).createChildPath("component_0");
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(componentPath));
  const auto& component0 = dataStructure.getDataRefAs<Float64Array>(componentPath);
  REQUIRE(component0.getNumberOfTuples() == k_ExpectedNumTuples);

  // The source .mat's ODFval is a normalised ODF (sums to 1.0) but is NOT a
  // per-cell uniform distribution — values vary across bins. Verify the two
  // properties we can actually assert: strict positivity of all bins and the
  // overall normalisation.
  double sum = 0.0;
  bool allPositive = true;
  for(usize i = 0; i < component0.getSize(); ++i)
  {
    const double v = component0[i];
    if(v <= 0.0)
    {
      allPositive = false;
    }
    sum += v;
  }
  REQUIRE(allPositive);
  REQUIRE(std::fabs(sum - 1.0) < 1e-9);
}

TEST_CASE("MTRSim::ReadMTRSimODFFilter: Rejects wrong path prefix", "[MTRSim][ReadMTRSimODFFilter][ErrorPath]")
{
  UnitTest::LoadPlugins();

  const fs::path inputFile = fs::path(fmt::format("{}/simulation_ODF.h5", unit_test::k_DataDir.view()));
  REQUIRE(fs::exists(inputFile));

  DataStructure dataStructure;

  ReadMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(inputFile));
  args.insertOrAssign(ReadMTRSimODFFilter::k_Hdf5PathPrefix_Key, std::make_any<StringParameter::ValueType>(std::string{"/does_not_exist"}));
  args.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key, std::make_any<DataObjectNameParameter::ValueType>(std::string{"Cell Data"}));

  auto preflightResult = filter.preflight(dataStructure, args);
  REQUIRE(preflightResult.outputActions.invalid());
}

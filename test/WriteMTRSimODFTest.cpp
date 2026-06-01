/**
 * Unit tests for WriteMTRSimODFFilter (Milestone AJ, Task 4 + Task 2
 * path-prefix cross-cut).
 *
 * Tests cover:
 *   1. Byte-exact round-trip through /ODF_best (back-compat). The prefix arg is
 *      passed explicitly to be sure the default still works end-to-end, and the
 *      new "HDF5 Path Prefix" PreflightUpdatedValue is asserted.
 *   2. An asymmetric-dims sentinel test that catches any axis-permutation bug
 *      in either Read or Write (closes the validation gap from the symmetric
 *      72x72 exemplar).
 *   3. Preflight error path: zero component arrays selected.
 *   4. Round-trip through a distinct prefix using the blank_ODF fixture.
 *   5. Round-trip through a distinct prefix using the uniform_ODF fixture.
 */

#include <catch2/catch.hpp>

#include "MTRSim/Filters/ReadMTRSimODFFilter.hpp"
#include "MTRSim/Filters/WriteMTRSimODFFilter.hpp"
#include "MTRSim/MTRSim_test_dirs.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/MultiArraySelectionParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include "LibMTRSim/ODFFileIO.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <system_error>

#if defined(_WIN32)
#include <process.h>
#define MTRSIM_GET_PID() static_cast<long long>(::_getpid())
#else
#include <unistd.h>
#define MTRSIM_GET_PID() static_cast<long long>(::getpid())
#endif

namespace fs = std::filesystem;

using namespace nx::core;
using namespace nx::core::UnitTest;

TEST_CASE(
    "MTRSim::WriteMTRSimODFFilter: Round-trip through /ODF_best (back-compat)",
    "[MTRSim][WriteMTRSimODFFilter]") {
  UnitTest::LoadPlugins();

  const fs::path inputFile = fs::path(
      fmt::format("{}/simulation_ODF.h5", unit_test::k_DataDir.view()));
  REQUIRE(fs::exists(inputFile));

  const fs::path outputFile =
      fs::temp_directory_path() /
      fmt::format("mtrsim_export_roundtrip_{}.h5", MTRSIM_GET_PID());
  std::error_code removeErr;
  if (fs::exists(outputFile)) {
    fs::remove(outputFile, removeErr);
  }

  const DataPath imageGeomPath({"ODF"});
  const std::string cellAttrMatName = "Cell Data";

  DataStructure dataStructure;

  // Import the reference file.
  {
    ReadMTRSimODFFilter readFilter;
    Arguments readArgs;
    readArgs.insertOrAssign(
        ReadMTRSimODFFilter::k_InputFile_Key,
        std::make_any<FileSystemPathParameter::ValueType>(inputFile));
    readArgs.insertOrAssign(
        ReadMTRSimODFFilter::k_Hdf5PathPrefix_Key,
        std::make_any<StringParameter::ValueType>(std::string{"/ODF_best"}));
    readArgs.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key,
                            std::make_any<DataPath>(imageGeomPath));
    readArgs.insertOrAssign(
        ReadMTRSimODFFilter::k_CellAttrMatName_Key,
        std::make_any<DataObjectNameParameter::ValueType>(cellAttrMatName));

    auto readPreflight = readFilter.preflight(dataStructure, readArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(readPreflight.outputActions);
    auto readExecute = readFilter.execute(dataStructure, readArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(readExecute.result);
  }

  // Determine the component paths that were created so we can pass them to the
  // writer.
  const auto origComponents = mtrsim::readODFComponents(inputFile);
  REQUIRE(origComponents.size() >= 1);

  const DataPath cellAttrMatPath =
      imageGeomPath.createChildPath(cellAttrMatName);
  std::vector<DataPath> componentPaths;
  componentPaths.reserve(origComponents.size());
  for (std::size_t c = 0; c < origComponents.size(); ++c) {
    componentPaths.push_back(
        cellAttrMatPath.createChildPath(fmt::format("component_{}", c)));
  }

  // Export to the temp file.
  {
    WriteMTRSimODFFilter writeFilter;
    Arguments writeArgs;
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_OutputFile_Key,
        std::make_any<FileSystemPathParameter::ValueType>(outputFile));
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_Hdf5PathPrefix_Key,
        std::make_any<StringParameter::ValueType>(std::string{"/ODF_best"}));
    writeArgs.insertOrAssign(WriteMTRSimODFFilter::k_InputImageGeometry_Key,
                             std::make_any<DataPath>(imageGeomPath));
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_ODFComponents_Key,
        std::make_any<MultiArraySelectionParameter::ValueType>(componentPaths));

    auto writePreflight = writeFilter.preflight(dataStructure, writeArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(writePreflight.outputActions);

    // Assert the prefix-preview updated value surfaced.
    const bool foundPrefixLabel = std::any_of(
        writePreflight.outputValues.begin(), writePreflight.outputValues.end(),
        [](const IFilter::PreflightValue &v) {
          return v.name == "HDF5 Path Prefix";
        });
    REQUIRE(foundPrefixLabel);

    auto writeExecute = writeFilter.execute(dataStructure, writeArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(writeExecute.result);
  }

  REQUIRE(fs::exists(outputFile));

  // Byte-exact comparison via libmtrsim.
  auto roundComponents = mtrsim::readODFComponents(outputFile);
  REQUIRE(origComponents.size() == roundComponents.size());
  for (std::size_t c = 0; c < origComponents.size(); ++c) {
    REQUIRE(origComponents[c].values == roundComponents[c].values);
  }

  fs::remove(outputFile, removeErr);
}

TEST_CASE(
    "MTRSim::WriteMTRSimODFFilter: axis mapping preserves phi1-PHI-phi2 layout",
    "[MTRSim][WriteMTRSimODFFilter][AxisMapping]") {
  UnitTest::LoadPlugins();

  // Asymmetric dimensions so an axis permutation cannot accidentally produce
  // the correct answer.
  constexpr std::size_t k_NumPhi1 = 24; // slowest on disk -> Z in ImageGeom
  constexpr std::size_t k_NumPHI = 18;  //                   -> Y in ImageGeom
  constexpr std::size_t k_NumPhi2 = 72; // fastest on disk -> X in ImageGeom
  constexpr float k_SpacingPhi1Deg = 15.0f;
  constexpr float k_SpacingPHIDeg = 10.0f;
  constexpr float k_SpacingPhi2Deg = 5.0f;

  const DataPath imageGeomPath({"ODF"});
  const std::string cellAttrMatName = "Cell Data";
  const DataPath cellAttrMatPath =
      imageGeomPath.createChildPath(cellAttrMatName);
  const DataPath componentPath = cellAttrMatPath.createChildPath("component_0");

  DataStructure dataStructure;

  // Build the ImageGeom via CreateImageGeometryAction (XYZ = phi2, PHI, phi1).
  {
    CreateImageGeometryAction::DimensionType dims = {k_NumPhi2, k_NumPHI,
                                                     k_NumPhi1};
    CreateImageGeometryAction::OriginType origin = {0.0f, 0.0f, 0.0f};
    CreateImageGeometryAction::SpacingType spacing = {
        k_SpacingPhi2Deg, k_SpacingPHIDeg, k_SpacingPhi1Deg};
    CreateImageGeometryAction geomAction(imageGeomPath, dims, origin, spacing,
                                         cellAttrMatName,
                                         IGeometry::LengthUnit::Unknown);
    auto geomResult =
        geomAction.apply(dataStructure, IDataAction::Mode::Execute);
    SIMPLNX_RESULT_REQUIRE_VALID(geomResult);
  }

  // Create the single Float64 component array with ZYX tuple shape {nPhi1,
  // nPHI, nPhi2}.
  {
    std::vector<usize> tupleShapeZYX = {k_NumPhi1, k_NumPHI, k_NumPhi2};
    CreateArrayAction arrayAction(DataType::float64, tupleShapeZYX,
                                  std::vector<usize>{1}, componentPath);
    auto arrayResult =
        arrayAction.apply(dataStructure, IDataAction::Mode::Execute);
    SIMPLNX_RESULT_REQUIRE_VALID(arrayResult);
  }

  // Initialize to zero, then place a single sentinel at a known (iPhi1, iPHI,
  // iPhi2).
  constexpr std::size_t k_iPhi1 = 3;
  constexpr std::size_t k_iPHI = 5;
  constexpr std::size_t k_iPhi2 = 7;
  // Row-major flat index with phi1 slowest, phi2 fastest: 3*18*72 + 5*72 + 7 =
  // 4255.
  constexpr std::size_t k_SentinelIndex =
      k_iPhi1 * (k_NumPHI * k_NumPhi2) + k_iPHI * k_NumPhi2 + k_iPhi2;
  static_assert(k_SentinelIndex == 4255,
                "Recompute sentinel flat index before asserting on it.");
  constexpr double k_SentinelValue = 42.0;

  {
    auto &componentArray =
        dataStructure.getDataRefAs<Float64Array>(componentPath);
    auto &store = componentArray.getDataStoreRef();
    for (std::size_t i = 0; i < store.getSize(); ++i) {
      store.setValue(i, 0.0);
    }
    store.setValue(k_SentinelIndex, k_SentinelValue);
  }

  // Write out via the filter.
  const fs::path outputFile =
      fs::temp_directory_path() /
      fmt::format("mtrsim_axis_mapping_{}.h5", MTRSIM_GET_PID());
  std::error_code removeErr;
  if (fs::exists(outputFile)) {
    fs::remove(outputFile, removeErr);
  }

  {
    WriteMTRSimODFFilter writeFilter;
    Arguments writeArgs;
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_OutputFile_Key,
        std::make_any<FileSystemPathParameter::ValueType>(outputFile));
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_Hdf5PathPrefix_Key,
        std::make_any<StringParameter::ValueType>(std::string{"/ODF_best"}));
    writeArgs.insertOrAssign(WriteMTRSimODFFilter::k_InputImageGeometry_Key,
                             std::make_any<DataPath>(imageGeomPath));
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_ODFComponents_Key,
        std::make_any<MultiArraySelectionParameter::ValueType>(
            std::vector<DataPath>{componentPath}));

    auto writePreflight = writeFilter.preflight(dataStructure, writeArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(writePreflight.outputActions);
    auto writeExecute = writeFilter.execute(dataStructure, writeArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(writeExecute.result);
  }

  REQUIRE(fs::exists(outputFile));

  // Read back the metadata and verify axis mapping was reversed correctly.
  auto meta = mtrsim::readODFMetadata(outputFile);
  REQUIRE(meta.numComponents == 1);
  REQUIRE(meta.dimsPhi1PHIPhi2 ==
          std::array<int64_t, 3>{static_cast<int64_t>(k_NumPhi1),
                                 static_cast<int64_t>(k_NumPHI),
                                 static_cast<int64_t>(k_NumPhi2)});
  REQUIRE(std::abs(meta.spacingDegPhi1PHIPhi2[0] -
                   static_cast<double>(k_SpacingPhi1Deg)) < 1e-9);
  REQUIRE(std::abs(meta.spacingDegPhi1PHIPhi2[1] -
                   static_cast<double>(k_SpacingPHIDeg)) < 1e-9);
  REQUIRE(std::abs(meta.spacingDegPhi1PHIPhi2[2] -
                   static_cast<double>(k_SpacingPhi2Deg)) < 1e-9);

  auto comps = mtrsim::readODFComponents(outputFile);
  REQUIRE(comps.size() == 1);
  const auto &vals = comps[0].values;
  REQUIRE(vals.size() == k_NumPhi1 * k_NumPHI * k_NumPhi2);
  REQUIRE(vals[k_SentinelIndex] == k_SentinelValue);

  std::size_t nonzeroCount = 0;
  for (std::size_t i = 0; i < vals.size(); ++i) {
    if (vals[i] != 0.0) {
      ++nonzeroCount;
    }
  }
  REQUIRE(nonzeroCount == 1);

  fs::remove(outputFile, removeErr);
}

TEST_CASE("MTRSim::WriteMTRSimODFFilter: zero components rejects at preflight",
          "[MTRSim][WriteMTRSimODFFilter][ErrorPath]") {
  UnitTest::LoadPlugins();

  DataStructure dataStructure;

  WriteMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(
      WriteMTRSimODFFilter::k_OutputFile_Key,
      std::make_any<FileSystemPathParameter::ValueType>(
          fs::temp_directory_path() / "mtrsim_should_not_be_written.h5"));
  args.insertOrAssign(
      WriteMTRSimODFFilter::k_Hdf5PathPrefix_Key,
      std::make_any<StringParameter::ValueType>(std::string{"/ODF_best"}));
  args.insertOrAssign(WriteMTRSimODFFilter::k_InputImageGeometry_Key,
                      std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(WriteMTRSimODFFilter::k_ODFComponents_Key,
                      std::make_any<MultiArraySelectionParameter::ValueType>(
                          std::vector<DataPath>{}));

  auto preflightResult = filter.preflight(dataStructure, args);
  REQUIRE(preflightResult.outputActions.invalid());
}

namespace {
// Helper: read fixture with a given prefix, write it out with a distinct
// prefix, read back and verify byte-exact equality of the component values.
void runDistinctPrefixRoundTrip(const std::string &fixtureName,
                                const std::string &readPrefix,
                                const std::string &writePrefix,
                                const std::string &tempStem) {
  UnitTest::LoadPlugins();

  const fs::path inputFile =
      fs::path(fmt::format("{}/{}", unit_test::k_DataDir.view(), fixtureName));
  REQUIRE(fs::exists(inputFile));

  const fs::path outputFile =
      fs::temp_directory_path() /
      fmt::format("{}_{}.h5", tempStem, MTRSIM_GET_PID());
  std::error_code removeErr;
  if (fs::exists(outputFile)) {
    fs::remove(outputFile, removeErr);
  }

  const DataPath imageGeomPath({"ODF"});
  const std::string cellAttrMatName = "Cell Data";

  DataStructure dataStructure;

  {
    ReadMTRSimODFFilter readFilter;
    Arguments readArgs;
    readArgs.insertOrAssign(
        ReadMTRSimODFFilter::k_InputFile_Key,
        std::make_any<FileSystemPathParameter::ValueType>(inputFile));
    readArgs.insertOrAssign(
        ReadMTRSimODFFilter::k_Hdf5PathPrefix_Key,
        std::make_any<StringParameter::ValueType>(readPrefix));
    readArgs.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key,
                            std::make_any<DataPath>(imageGeomPath));
    readArgs.insertOrAssign(
        ReadMTRSimODFFilter::k_CellAttrMatName_Key,
        std::make_any<DataObjectNameParameter::ValueType>(cellAttrMatName));

    auto readPreflight = readFilter.preflight(dataStructure, readArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(readPreflight.outputActions);
    auto readExecute = readFilter.execute(dataStructure, readArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(readExecute.result);
  }

  const auto origComponents = mtrsim::readODFComponents(inputFile, readPrefix);
  REQUIRE(origComponents.size() >= 1);

  const DataPath cellAttrMatPath =
      imageGeomPath.createChildPath(cellAttrMatName);
  std::vector<DataPath> componentPaths;
  componentPaths.reserve(origComponents.size());
  for (std::size_t c = 0; c < origComponents.size(); ++c) {
    componentPaths.push_back(
        cellAttrMatPath.createChildPath(fmt::format("component_{}", c)));
  }

  {
    WriteMTRSimODFFilter writeFilter;
    Arguments writeArgs;
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_OutputFile_Key,
        std::make_any<FileSystemPathParameter::ValueType>(outputFile));
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_Hdf5PathPrefix_Key,
        std::make_any<StringParameter::ValueType>(writePrefix));
    writeArgs.insertOrAssign(WriteMTRSimODFFilter::k_InputImageGeometry_Key,
                             std::make_any<DataPath>(imageGeomPath));
    writeArgs.insertOrAssign(
        WriteMTRSimODFFilter::k_ODFComponents_Key,
        std::make_any<MultiArraySelectionParameter::ValueType>(componentPaths));

    auto writePreflight = writeFilter.preflight(dataStructure, writeArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(writePreflight.outputActions);
    auto writeExecute = writeFilter.execute(dataStructure, writeArgs);
    SIMPLNX_RESULT_REQUIRE_VALID(writeExecute.result);
  }

  REQUIRE(fs::exists(outputFile));

  // The default prefix should NOT exist at the destination — only writePrefix.
  // Library reads against the default "/ODF_best" must fail.
  // (Only meaningful when writePrefix != "/ODF_best".)
  if (writePrefix != "/ODF_best") {
    REQUIRE_THROWS(mtrsim::readODFMetadata(outputFile));
  }

  auto roundComponents = mtrsim::readODFComponents(outputFile, writePrefix);
  REQUIRE(origComponents.size() == roundComponents.size());
  for (std::size_t c = 0; c < origComponents.size(); ++c) {
    REQUIRE(origComponents[c].values == roundComponents[c].values);
  }

  fs::remove(outputFile, removeErr);
}
} // namespace

TEST_CASE(
    "MTRSim::WriteMTRSimODFFilter: Round-trip blank_ODF with distinct prefix",
    "[MTRSim][WriteMTRSimODFFilter]") {
  runDistinctPrefixRoundTrip("blank_ODF.h5", "/blank_ODF",
                             "/blank_ODF_exported", "mtrsim_blank_rt");
}

TEST_CASE(
    "MTRSim::WriteMTRSimODFFilter: Round-trip uniform_ODF with distinct prefix",
    "[MTRSim][WriteMTRSimODFFilter]") {
  runDistinctPrefixRoundTrip("uniform_ODF.h5", "/uniform_ODF",
                             "/uniform_ODF_exported", "mtrsim_uniform_rt");
}

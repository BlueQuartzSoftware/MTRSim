#include <catch2/catch.hpp>

#include "MTRSim/Filters/WriteODFSectionsFilter.hpp"

#include "EbsdLib/Core/EbsdLibConstants.h"
#include "EbsdLib/Utilities/ODFSectionCompositor.h"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataStore.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/StringArray.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/CreateColorMapParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataPathSelectionParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <limits>
#include <string>
#include <vector>

using namespace nx::core;

namespace
{
constexpr usize k_NumPhi1 = 72;
constexpr usize k_NumPhi = 36;
constexpr usize k_NumPhi2 = 72;

const DataPath k_OdfGeometryPath({"ODF"});
const DataPath k_OdfArrayPath({"ODF", "Cell Data", "ODF Values"});
const DataPath k_CrystalStructuresPath({"Ensemble Data", "Crystal Structures"});
const DataPath k_MaterialNamesPath({"Ensemble Data", "Material Names"});
const DataPath k_OutputGeometryPath({"ODF Sections"});

DataStructure CreateValidDataStructure()
{
  DataStructure dataStructure;

  auto* imageGeom = ImageGeom::Create(dataStructure, k_OdfGeometryPath.getTargetName());
  imageGeom->setDimensions({k_NumPhi2, k_NumPhi, k_NumPhi1});
  imageGeom->setOrigin({0.0F, 0.0F, 0.0F});
  imageGeom->setSpacing({5.0F, 5.0F, 5.0F});

  const ShapeType tupleShape = {k_NumPhi1, k_NumPhi, k_NumPhi2};
  auto* cellData = AttributeMatrix::Create(dataStructure, "Cell Data", tupleShape, imageGeom->getId());
  imageGeom->setCellData(cellData->getId());
  Float64Array::CreateWithStore<DataStore<float64>>(dataStructure, k_OdfArrayPath.getTargetName(), tupleShape, {1}, cellData->getId());

  auto* ensembleData = AttributeMatrix::Create(dataStructure, "Ensemble Data", {2});
  auto* crystalStructures = UInt32Array::CreateWithStore<DataStore<uint32>>(dataStructure, k_CrystalStructuresPath.getTargetName(), {2}, {1}, ensembleData->getId());
  (*crystalStructures)[0] = ebsdlib::CrystalStructure::UnknownCrystalStructure;
  (*crystalStructures)[1] = ebsdlib::CrystalStructure::Hexagonal_High;
  StringArray::CreateWithValues(dataStructure, k_MaterialNamesPath.getTargetName(), {2}, {"Invalid", "Titanium (Alpha)"}, ensembleData->getId());
  return dataStructure;
}

Arguments CreateValidArguments()
{
  Arguments args;
  args.insertOrAssign(WriteODFSectionsFilter::k_InputImageGeometryPath_Key, k_OdfGeometryPath);
  args.insertOrAssign(WriteODFSectionsFilter::k_ODFArrayPath_Key, k_OdfArrayPath);
  args.insertOrAssign(WriteODFSectionsFilter::k_InputUnits_Key, ChoicesParameter::ValueType{0});
  args.insertOrAssign(WriteODFSectionsFilter::k_CrystalStructuresPath_Key, k_CrystalStructuresPath);
  args.insertOrAssign(WriteODFSectionsFilter::k_MaterialNamesPath_Key, k_MaterialNamesPath);
  args.insertOrAssign(WriteODFSectionsFilter::k_PhaseIndex_Key, int32{1});
  args.insertOrAssign(WriteODFSectionsFilter::k_Title_Key, std::string{"ODF Sections"});
  args.insertOrAssign(WriteODFSectionsFilter::k_NumSections_Key, int32{6});
  args.insertOrAssign(WriteODFSectionsFilter::k_SectionsPerRow_Key, int32{3});
  args.insertOrAssign(WriteODFSectionsFilter::k_SectionWidth_Key, int32{512});
  args.insertOrAssign(WriteODFSectionsFilter::k_ColorTablePreset_Key, CreateColorMapParameter::ValueType{"Rainbow"});
  args.insertOrAssign(WriteODFSectionsFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{0});
  args.insertOrAssign(WriteODFSectionsFilter::k_ManualMinimum_Key, float64{0.0});
  args.insertOrAssign(WriteODFSectionsFilter::k_ManualMaximum_Key, float64{1.0});
  args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, true);
  args.insertOrAssign(WriteODFSectionsFilter::k_OutputPngPath_Key, FileSystemPathParameter::ValueType{"odf_sections.png"});
  args.insertOrAssign(WriteODFSectionsFilter::k_SaveAsImageGeometry_Key, true);
  args.insertOrAssign(WriteODFSectionsFilter::k_OutputImageGeometryPath_Key, k_OutputGeometryPath);
  return args;
}

std::string GetPreflightValue(const IFilter::PreflightResult& result, const std::string& name)
{
  const auto iter = std::find_if(result.outputValues.begin(), result.outputValues.end(), [&name](const IFilter::PreflightValue& value) { return value.name == name; });
  REQUIRE(iter != result.outputValues.end());
  return iter->value;
}

struct InvalidCase
{
  std::string name;
  int32 errorCode = 0;
  std::string diagnosticFragment;
  std::function<void(DataStructure&, Arguments&)> mutate;
};
} // namespace

TEST_CASE("MTRSim::WriteODFSectionsFilter: Parameter Preflight exposes the standard ODF section contract", "[MTRSim][WriteODFSectionsFilter][Preflight][Parameters]")
{
  WriteODFSectionsFilter filter;
  REQUIRE(filter.uuid() == *Uuid::FromString("5ec65380-b1a8-4aff-8bcd-310b5067b236"));

  const Parameters params = filter.parameters();
  REQUIRE(params.size() == 18);
  REQUIRE(dynamic_cast<const GeometrySelectionParameter*>(params.at(WriteODFSectionsFilter::k_InputImageGeometryPath_Key).get()) != nullptr);
  const auto* odfArrayParameter = dynamic_cast<const ArraySelectionParameter*>(params.at(WriteODFSectionsFilter::k_ODFArrayPath_Key).get());
  REQUIRE(odfArrayParameter != nullptr);
  const ArraySelectionParameter::AllowedTypes expectedOdfTypes = {DataType::float32, DataType::float64};
  const ArraySelectionParameter::AllowedComponentShapes expectedOdfComponentShapes = {{1}};
  REQUIRE(odfArrayParameter->allowedTypes() == expectedOdfTypes);
  REQUIRE(odfArrayParameter->requiredComponentShapes() == expectedOdfComponentShapes);
  REQUIRE(dynamic_cast<const ChoicesParameter*>(params.at(WriteODFSectionsFilter::k_InputUnits_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const ArraySelectionParameter*>(params.at(WriteODFSectionsFilter::k_CrystalStructuresPath_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const DataPathSelectionParameter*>(params.at(WriteODFSectionsFilter::k_MaterialNamesPath_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const Int32Parameter*>(params.at(WriteODFSectionsFilter::k_PhaseIndex_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const StringParameter*>(params.at(WriteODFSectionsFilter::k_Title_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const Int32Parameter*>(params.at(WriteODFSectionsFilter::k_NumSections_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const Int32Parameter*>(params.at(WriteODFSectionsFilter::k_SectionsPerRow_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const Int32Parameter*>(params.at(WriteODFSectionsFilter::k_SectionWidth_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const CreateColorMapParameter*>(params.at(WriteODFSectionsFilter::k_ColorTablePreset_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const ChoicesParameter*>(params.at(WriteODFSectionsFilter::k_ScaleMode_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const Float64Parameter*>(params.at(WriteODFSectionsFilter::k_ManualMinimum_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const Float64Parameter*>(params.at(WriteODFSectionsFilter::k_ManualMaximum_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const BoolParameter*>(params.at(WriteODFSectionsFilter::k_WritePng_Key).get()) != nullptr);
  const auto* outputPngParameter = dynamic_cast<const FileSystemPathParameter*>(params.at(WriteODFSectionsFilter::k_OutputPngPath_Key).get());
  REQUIRE(outputPngParameter != nullptr);
  REQUIRE(outputPngParameter->getAvailableExtensions() == FileSystemPathParameter::ExtensionsType{".png"});
  REQUIRE(outputPngParameter->getPathType() == FileSystemPathParameter::PathType::OutputFile);
  REQUIRE_FALSE(outputPngParameter->acceptAllExtensions());
  REQUIRE(dynamic_cast<const BoolParameter*>(params.at(WriteODFSectionsFilter::k_SaveAsImageGeometry_Key).get()) != nullptr);
  REQUIRE(dynamic_cast<const DataGroupCreationParameter*>(params.at(WriteODFSectionsFilter::k_OutputImageGeometryPath_Key).get()) != nullptr);

  const Arguments defaults = filter.getDefaultArguments();
  REQUIRE(defaults.value<ChoicesParameter::ValueType>(WriteODFSectionsFilter::k_InputUnits_Key) == 0);
  REQUIRE(defaults.value<int32>(WriteODFSectionsFilter::k_PhaseIndex_Key) == 1);
  REQUIRE(defaults.value<int32>(WriteODFSectionsFilter::k_NumSections_Key) == 6);
  REQUIRE(defaults.value<int32>(WriteODFSectionsFilter::k_SectionsPerRow_Key) == 3);
  REQUIRE(defaults.value<int32>(WriteODFSectionsFilter::k_SectionWidth_Key) == 512);
  REQUIRE(defaults.value<CreateColorMapParameter::ValueType>(WriteODFSectionsFilter::k_ColorTablePreset_Key) == "Rainbow");
  REQUIRE(defaults.value<ChoicesParameter::ValueType>(WriteODFSectionsFilter::k_ScaleMode_Key) == 0);
  REQUIRE(defaults.value<float64>(WriteODFSectionsFilter::k_ManualMinimum_Key) == 0.0);
  REQUIRE(defaults.value<float64>(WriteODFSectionsFilter::k_ManualMaximum_Key) == 1.0);
  REQUIRE(defaults.value<bool>(WriteODFSectionsFilter::k_WritePng_Key));
  REQUIRE(defaults.value<bool>(WriteODFSectionsFilter::k_SaveAsImageGeometry_Key));
  REQUIRE(defaults.value<DataPath>(WriteODFSectionsFilter::k_OutputImageGeometryPath_Key) == k_OutputGeometryPath);
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Valid Preflight creates the composite image", "[MTRSim][WriteODFSectionsFilter][Preflight]")
{
  DataStructure dataStructure = CreateValidDataStructure();
  WriteODFSectionsFilter filter;
  const auto result = filter.preflight(dataStructure, CreateValidArguments());
  SIMPLNX_RESULT_REQUIRE_VALID(result.outputActions);

  REQUIRE(GetPreflightValue(result, "Section Angles") == "0, 10, 20, 30, 40, 50 degrees");
  REQUIRE(GetPreflightValue(result, "Phase") == "1");
  REQUIRE(GetPreflightValue(result, "Material") == "Titanium (Alpha)");
  REQUIRE(GetPreflightValue(result, "Laue Class") == "Hexagonal 6/mmm (D6h)");
  REQUIRE(GetPreflightValue(result, "Input Units") == "MUD");
  REQUIRE(GetPreflightValue(result, "Rendered Units") == "MUD");
  REQUIRE(GetPreflightValue(result, "Full Grid Dimensions (phi1 x PHI x phi2)") == "72 x 36 x 72");
  REQUIRE(GetPreflightValue(result, "Full Grid Spacing (phi1, PHI, phi2)") == "5, 5, 5 degrees");
  REQUIRE(GetPreflightValue(result, "Section Layout") == "3 columns x 2 rows");
  REQUIRE(GetPreflightValue(result, "Scale Mode") == "Automatic");

  const auto& actions = result.outputActions.value().actions;
  REQUIRE(actions.size() == 2);

  const auto* imageAction = dynamic_cast<const CreateImageGeometryAction*>(actions[0].get());
  REQUIRE(imageAction != nullptr);
  const auto* arrayAction = dynamic_cast<const CreateArrayAction*>(actions[1].get());
  REQUIRE(arrayAction != nullptr);

  ebsdlib::ODFSectionConfiguration config;
  config.grid.laueOpsIndex = ebsdlib::CrystalStructure::Hexagonal_High;
  config.sectionCount = 6;
  config.sectionsPerRow = 3;
  config.sectionWidth = 512;
  const auto layout = ebsdlib::ComputeODFSectionLayout(config);

  REQUIRE(imageAction->path() == k_OutputGeometryPath);
  REQUIRE(imageAction->dims() == CreateImageGeometryAction::DimensionType{static_cast<usize>(layout.pageWidth), static_cast<usize>(layout.pageHeight), 1});
  REQUIRE(imageAction->cellAttributeMatrixName() == "Cell Data");

  REQUIRE(arrayAction->path() == DataPath({"ODF Sections", "Cell Data", "ODF Sections RGB"}));
  REQUIRE(arrayAction->type() == DataType::uint8);
  REQUIRE(arrayAction->dims() == ShapeType{1, static_cast<usize>(layout.pageHeight), static_cast<usize>(layout.pageWidth)});
  REQUIRE(arrayAction->componentDims() == ShapeType{3});
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: PNG-only Preflight creates no DataStructure actions", "[MTRSim][WriteODFSectionsFilter][Preflight]")
{
  DataStructure dataStructure = CreateValidDataStructure();
  Arguments args = CreateValidArguments();
  args.insertOrAssign(WriteODFSectionsFilter::k_SaveAsImageGeometry_Key, false);

  WriteODFSectionsFilter filter;
  const auto result = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(result.outputActions);
  REQUIRE(result.outputActions.value().actions.empty());
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Preflight previews both unit and scale modes", "[MTRSim][WriteODFSectionsFilter][Preflight]")
{
  struct PreviewCase
  {
    ChoicesParameter::ValueType inputUnits = 0;
    ChoicesParameter::ValueType scaleMode = 0;
    std::string expectedInputUnits;
    std::string expectedScaleMode;
  };

  const std::vector<PreviewCase> cases = {
      {0, 0, "MUD", "Automatic"},
      {1, 0, "Count-Density", "Automatic"},
      {0, 1, "MUD", "Manual: [0.25, 2.5] MUD"},
      {1, 1, "Count-Density", "Manual: [0.25, 2.5] MUD"},
  };

  for(const auto& testCase : cases)
  {
    DYNAMIC_SECTION("units " << testCase.expectedInputUnits << ", scale " << testCase.expectedScaleMode)
    {
      DataStructure dataStructure = CreateValidDataStructure();
      Arguments args = CreateValidArguments();
      args.insertOrAssign(WriteODFSectionsFilter::k_InputUnits_Key, testCase.inputUnits);
      args.insertOrAssign(WriteODFSectionsFilter::k_ScaleMode_Key, testCase.scaleMode);
      args.insertOrAssign(WriteODFSectionsFilter::k_ManualMinimum_Key, float64{0.25});
      args.insertOrAssign(WriteODFSectionsFilter::k_ManualMaximum_Key, float64{2.5});

      WriteODFSectionsFilter filter;
      const auto result = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(result.outputActions);
      REQUIRE(GetPreflightValue(result, "Input Units") == testCase.expectedInputUnits);
      REQUIRE(GetPreflightValue(result, "Rendered Units") == "MUD");
      REQUIRE(GetPreflightValue(result, "Scale Mode") == testCase.expectedScaleMode);
    }
  }
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Extreme section count Preflight rejects layout without allocating angles", "[MTRSim][WriteODFSectionsFilter][Preflight][ErrorPath]")
{
  DataStructure dataStructure = CreateValidDataStructure();
  Arguments args = CreateValidArguments();
  args.insertOrAssign(WriteODFSectionsFilter::k_NumSections_Key, std::numeric_limits<int32>::max());
  args.insertOrAssign(WriteODFSectionsFilter::k_SectionsPerRow_Key, int32{1});

  WriteODFSectionsFilter filter;
  const auto start = std::chrono::steady_clock::now();
  const auto result = filter.preflight(dataStructure, args);
  const auto elapsed = std::chrono::steady_clock::now() - start;
  SIMPLNX_RESULT_REQUIRE_INVALID(result.outputActions);
  REQUIRE(result.outputActions.errors()[0].code == -12315);
  REQUIRE(result.outputActions.errors()[0].message.find(std::to_string(std::numeric_limits<int32>::max())) != std::string::npos);
  REQUIRE(elapsed < std::chrono::seconds{1});
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Invalid Preflight reports actionable errors", "[MTRSim][WriteODFSectionsFilter][Preflight][ErrorPath]")
{
  const std::vector<InvalidCase> cases = {
      {"both outputs disabled", -12300, "false",
       [](DataStructure&, Arguments& args) {
         args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, false);
         args.insertOrAssign(WriteODFSectionsFilter::k_SaveAsImageGeometry_Key, false);
       }},
      {"ODF type", FilterParameter::Constants::k_Validate_AllowedType_Error, DataPath({"ODF", "Cell Data", "Integer ODF"}).toString(),
       [](DataStructure& dataStructure, Arguments& args) {
         const auto* cellData = dataStructure.getDataAs<AttributeMatrix>(DataPath({"ODF", "Cell Data"}));
         UInt32Array::CreateWithStore<DataStore<uint32>>(dataStructure, "Integer ODF", {k_NumPhi1, k_NumPhi, k_NumPhi2}, {1}, cellData->getId());
         args.insertOrAssign(WriteODFSectionsFilter::k_ODFArrayPath_Key, DataPath({"ODF", "Cell Data", "Integer ODF"}));
       }},
      {"ODF component shape", FilterParameter::Constants::k_Validate_TupleShapeValue, DataPath({"ODF", "Cell Data", "Vector ODF"}).toString(),
       [](DataStructure& dataStructure, Arguments& args) {
         const auto* cellData = dataStructure.getDataAs<AttributeMatrix>(DataPath({"ODF", "Cell Data"}));
         Float64Array::CreateWithStore<DataStore<float64>>(dataStructure, "Vector ODF", {k_NumPhi1, k_NumPhi, k_NumPhi2}, {2}, cellData->getId());
         args.insertOrAssign(WriteODFSectionsFilter::k_ODFArrayPath_Key, DataPath({"ODF", "Cell Data", "Vector ODF"}));
       }},
      {"ODF parent", -12302, DataPath({"Other Cell Data", "ODF Values"}).toString(),
       [](DataStructure& dataStructure, Arguments& args) {
         auto* otherData = AttributeMatrix::Create(dataStructure, "Other Cell Data", {k_NumPhi1, k_NumPhi, k_NumPhi2});
         Float64Array::CreateWithStore<DataStore<float64>>(dataStructure, "ODF Values", {k_NumPhi1, k_NumPhi, k_NumPhi2}, {1}, otherData->getId());
         args.insertOrAssign(WriteODFSectionsFilter::k_ODFArrayPath_Key, DataPath({"Other Cell Data", "ODF Values"}));
       }},
      {"nonzero origin", -12303, "1",
       [](DataStructure& dataStructure, Arguments&) {
         auto* imageGeom = dataStructure.getDataAs<ImageGeom>(k_OdfGeometryPath);
         imageGeom->setOrigin({0.0F, 1.0F, 0.0F});
       }},
      {"nonuniform spacing", -12304, "6",
       [](DataStructure& dataStructure, Arguments&) {
         auto* imageGeom = dataStructure.getDataAs<ImageGeom>(k_OdfGeometryPath);
         imageGeom->setSpacing({5.0F, 6.0F, 5.0F});
       }},
      {"zero spacing", -12304, "0",
       [](DataStructure& dataStructure, Arguments&) {
         auto* imageGeom = dataStructure.getDataAs<ImageGeom>(k_OdfGeometryPath);
         imageGeom->setSpacing({0.0F, 0.0F, 0.0F});
       }},
      {"negative spacing", -12304, "-5",
       [](DataStructure& dataStructure, Arguments&) {
         auto* imageGeom = dataStructure.getDataAs<ImageGeom>(k_OdfGeometryPath);
         imageGeom->setSpacing({-5.0F, -5.0F, -5.0F});
       }},
      {"Euler cube dimensions", -12305, "71",
       [](DataStructure& dataStructure, Arguments&) {
         auto* imageGeom = dataStructure.getDataAs<ImageGeom>(k_OdfGeometryPath);
         imageGeom->setDimensions({71, k_NumPhi, k_NumPhi1});
       }},
      {"ensemble tuple mismatch", -12306, "3",
       [](DataStructure& dataStructure, Arguments& args) {
         auto* ensembleData = AttributeMatrix::Create(dataStructure, "Three Ensemble Data", {3});
         StringArray::CreateWithValues(dataStructure, "Material Names", {3}, {"Invalid", "Titanium (Alpha)", "Extra"}, ensembleData->getId());
         args.insertOrAssign(WriteODFSectionsFilter::k_MaterialNamesPath_Key, DataPath({"Three Ensemble Data", "Material Names"}));
       }},
      {"phase zero", -12307, "0", [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_PhaseIndex_Key, int32{0}); }},
      {"phase out of range", -12307, "2", [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_PhaseIndex_Key, int32{2}); }},
      {"unknown crystal structure", -12308, "999",
       [](DataStructure& dataStructure, Arguments&) {
         auto* crystalStructures = dataStructure.getDataAs<UInt32Array>(k_CrystalStructuresPath);
         (*crystalStructures)[1] = ebsdlib::CrystalStructure::UnknownCrystalStructure;
       }},
      {"section count", -12309, "1", [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_NumSections_Key, int32{1}); }},
      {"sections per row zero", -12310, "0", [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_SectionsPerRow_Key, int32{0}); }},
      {"sections per row too large", -12310, "7", [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_SectionsPerRow_Key, int32{7}); }},
      {"section width", -12311, "127", [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_SectionWidth_Key, int32{127}); }},
      {"invalid preset", -12312, "Not A Preset",
       [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_ColorTablePreset_Key, CreateColorMapParameter::ValueType{"Not A Preset"}); }},
      {"negative manual minimum", -12313, "-0.5",
       [](DataStructure&, Arguments& args) {
         args.insertOrAssign(WriteODFSectionsFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
         args.insertOrAssign(WriteODFSectionsFilter::k_ManualMinimum_Key, float64{-0.5});
       }},
      {"non-finite manual minimum", -12313, "nan",
       [](DataStructure&, Arguments& args) {
         args.insertOrAssign(WriteODFSectionsFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
         args.insertOrAssign(WriteODFSectionsFilter::k_ManualMinimum_Key, std::numeric_limits<float64>::quiet_NaN());
       }},
      {"equal manual limits", -12313, "[1, 1]",
       [](DataStructure&, Arguments& args) {
         args.insertOrAssign(WriteODFSectionsFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
         args.insertOrAssign(WriteODFSectionsFilter::k_ManualMinimum_Key, float64{1.0});
         args.insertOrAssign(WriteODFSectionsFilter::k_ManualMaximum_Key, float64{1.0});
       }},
      {"reversed manual limits", -12313, "[2, 1]",
       [](DataStructure&, Arguments& args) {
         args.insertOrAssign(WriteODFSectionsFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
         args.insertOrAssign(WriteODFSectionsFilter::k_ManualMinimum_Key, float64{2.0});
         args.insertOrAssign(WriteODFSectionsFilter::k_ManualMaximum_Key, float64{1.0});
       }},
      {"non-PNG extension", -3003, ".jpg",
       [](DataStructure&, Arguments& args) { args.insertOrAssign(WriteODFSectionsFilter::k_OutputPngPath_Key, FileSystemPathParameter::ValueType{"odf_sections.jpg"}); }},
  };

  for(const auto& testCase : cases)
  {
    DYNAMIC_SECTION(testCase.name)
    {
      DataStructure dataStructure = CreateValidDataStructure();
      Arguments args = CreateValidArguments();
      testCase.mutate(dataStructure, args);

      WriteODFSectionsFilter filter;
      const auto result = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_INVALID(result.outputActions);
      REQUIRE_FALSE(result.outputActions.errors().empty());
      CHECK(result.outputActions.errors()[0].code == testCase.errorCode);
      CHECK(result.outputActions.errors()[0].message.find(testCase.diagnosticFragment) != std::string::npos);
    }
  }
}

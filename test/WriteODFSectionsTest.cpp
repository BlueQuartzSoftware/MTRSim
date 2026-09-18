#include <catch2/catch.hpp>

#include "MTRSim/Filters/Algorithms/WriteODFSections.hpp"
#include "MTRSim/Filters/ComputeODFFilter.hpp"
#include "MTRSim/Filters/WriteODFSectionsFilter.hpp"
#include "SimplnxCore/Filters/ReadImageFilter.hpp"

#include "EbsdLib/Core/EbsdLibConstants.h"
#include "EbsdLib/Utilities/ODFSectionChrome.h"
#include "EbsdLib/Utilities/ODFSectionCompositor.h"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataStore.hpp"
#include "simplnx/DataStructure/EmptyDataStore.hpp"
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
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <numbers>
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

struct TemporaryOutput
{
  const std::filesystem::path directory = std::filesystem::temp_directory_path() / ("odf-sections-" + Uuid::GenerateV4().str());
  const std::filesystem::path png = directory / "sections.png";

  TemporaryOutput()
  {
    REQUIRE(std::filesystem::create_directory(directory));
  }

  ~TemporaryOutput()
  {
    std::error_code error;
    std::filesystem::remove_all(directory, error);
  }
};

const DataPath k_RgbPath = k_OutputGeometryPath.createChildPath("Cell Data").createChildPath("ODF Sections RGB");

Arguments ExecutionArguments()
{
  auto args = CreateValidArguments();
  args.insertOrAssign(WriteODFSectionsFilter::k_SectionWidth_Key, int32{128});
  args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, false);
  return args;
}

void FillOdf(DataStructure& dataStructure, double value)
{
  auto* array = dataStructure.getDataAs<Float64Array>(k_OdfArrayPath);
  REQUIRE(array != nullptr);
  std::fill(array->begin(), array->end(), value);
}

DataStructure ReadPng(const std::filesystem::path& path)
{
  ReadImageFilter filter;
  auto args = filter.getDefaultArguments();
  args.insertOrAssign(ReadImageFilter::k_FileName_Key, path);
  args.insertOrAssign(ReadImageFilter::k_ImageGeometryPath_Key, DataPath({"Decoded"}));
  args.insertOrAssign(ReadImageFilter::k_CellDataName_Key, std::string{"Cells"});
  args.insertOrAssign(ReadImageFilter::k_ImageDataArrayPath_Key, std::string{"Pixels"});
  DataStructure data;
  const auto result = filter.execute(data, args);
  SIMPLNX_RESULT_REQUIRE_VALID(result.result);
  return data;
}

DataStructure CreateDeferredMetadataDataStructure()
{
  auto data = CreateValidDataStructure();
  REQUIRE(data.removeData(k_CrystalStructuresPath));
  REQUIRE(data.removeData(k_MaterialNamesPath));
  const CreateArrayAction action(DataType::uint32, {2}, {1}, k_CrystalStructuresPath);
  const auto result = action.apply(data, IDataAction::Mode::Preflight);
  SIMPLNX_RESULT_REQUIRE_VALID(result);
  const auto* ensemble = data.getDataAs<AttributeMatrix>(k_CrystalStructuresPath.getParent());
  REQUIRE(ensemble != nullptr);
  REQUIRE(StringArray::CreateWithValues(data, k_MaterialNamesPath.getTargetName(), {2}, {"", ""}, ensemble->getId()) != nullptr);
  return data;
}

WriteODFSectionsInputValues ExecutionInputValues()
{
  WriteODFSectionsInputValues values;
  values.inputImageGeometryPath = k_OdfGeometryPath;
  values.odfArrayPath = k_OdfArrayPath;
  values.crystalStructuresPath = k_CrystalStructuresPath;
  values.materialNamesPath = k_MaterialNamesPath;
  values.sectionWidth = 128;
  values.writePng = false;
  return values;
}
} // namespace

TEST_CASE("MTRSim::WriteODFSectionsFilter: Parameter Preflight exposes the standard ODF section contract", "[MTRSim][WriteODFSectionsFilter][Preflight][Parameters]")
{
  WriteODFSectionsFilter filter;
  REQUIRE(filter.uuid() == *Uuid::FromString("5ec65380-b1a8-4aff-8bcd-310b5067b236"));
  REQUIRE(filter.humanName() == "Write ODF Section Image");

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

TEST_CASE("MTRSim::WriteODFSectionsFilter: Float32 fractional spacing agrees with full grid validation", "[MTRSim][WriteODFSectionsFilter][Preflight]")
{
  const usize phiCount = GENERATE(usize{7}, usize{53});
  const float32 spacing = static_cast<float32>(180.0 / static_cast<double>(phiCount));
  auto data = CreateValidDataStructure();
  auto* geometry = data.getDataAs<ImageGeom>(k_OdfGeometryPath);
  REQUIRE(geometry != nullptr);
  geometry->setDimensions({2 * phiCount, phiCount, 2 * phiCount});
  geometry->setSpacing({spacing, spacing, spacing});
  geometry->getCellData()->resizeTuples({2 * phiCount, phiCount, 2 * phiCount});
  WriteODFSectionsFilter filter;
  const auto valid = filter.preflight(data, ExecutionArguments());
  SIMPLNX_RESULT_REQUIRE_VALID(valid.outputActions);
  const float32 invalidSpacing = spacing * 1.000001F;
  geometry->setSpacing({invalidSpacing, invalidSpacing, invalidSpacing});
  const auto invalid = filter.preflight(data, ExecutionArguments());
  REQUIRE(invalid.outputActions.invalid());
  REQUIRE(invalid.outputActions.errors().front().code == -12305);
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

TEST_CASE("MTRSim::WriteODFSectionsFilter: Units produce identical MUD and Count-Density section pixels", "[MTRSim][WriteODFSectionsFilter][Units]")
{
  auto mudData = CreateValidDataStructure();
  auto countData = CreateValidDataStructure();
  auto* mud = mudData.getDataAs<Float64Array>(k_OdfArrayPath);
  auto* countDensity = countData.getDataAs<Float64Array>(k_OdfArrayPath);
  REQUIRE(mud != nullptr);
  REQUIRE(countDensity != nullptr);
  const double stepRad = 5.0 * std::numbers::pi / 180.0;
  for(usize index = 0; index < mud->getSize(); ++index)
  {
    const usize phiIndex = (index / k_NumPhi2) % k_NumPhi;
    const double phiCenter = (static_cast<double>(phiIndex) + 0.5) * stepRad;
    (*mud)[index] = 1.0 + static_cast<double>((index / k_NumPhi2) % 7) * 0.125;
    (*countDensity)[index] = (*mud)[index] * stepRad * stepRad * stepRad * std::sin(phiCenter) / (8.0 * std::numbers::pi * std::numbers::pi);
  }
  auto args = CreateValidArguments();
  args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, false);
  args.insertOrAssign(WriteODFSectionsFilter::k_SectionWidth_Key, int32{128});
  WriteODFSectionsFilter filter;
  const auto mudResult = filter.execute(mudData, args);
  SIMPLNX_RESULT_REQUIRE_VALID(mudResult.result);
  args.insertOrAssign(WriteODFSectionsFilter::k_InputUnits_Key, ChoicesParameter::ValueType{1});
  const auto countResult = filter.execute(countData, args);
  SIMPLNX_RESULT_REQUIRE_VALID(countResult.result);
  const auto* mudGeom = mudData.getDataAs<ImageGeom>(k_OutputGeometryPath);
  const auto* countGeom = countData.getDataAs<ImageGeom>(k_OutputGeometryPath);
  REQUIRE(mudGeom != nullptr);
  REQUIRE(countGeom != nullptr);
  REQUIRE(mudGeom->getDimensions() == countGeom->getDimensions());
  const auto rgbPath = k_OutputGeometryPath.createChildPath("Cell Data").createChildPath("ODF Sections RGB");
  const auto* mudRgb = mudData.getDataAs<UInt8Array>(rgbPath);
  const auto* countRgb = countData.getDataAs<UInt8Array>(rgbPath);
  REQUIRE(mudRgb != nullptr);
  REQUIRE(countRgb != nullptr);
  ebsdlib::ODFSectionConfiguration config;
  config.grid.laueOpsIndex = ebsdlib::CrystalStructure::Hexagonal_High;
  config.sectionWidth = 128;
  const auto layout = ebsdlib::ComputeODFSectionLayout(config);
  const auto requireEqualRegion = [&](usize x, usize y, usize width, usize height) {
    bool equal = true;
    for(usize row = y; row < y + height; ++row)
    {
      for(usize column = x; column < x + width; ++column)
      {
        for(usize compIdx = 0; compIdx < 3; ++compIdx)
        {
          const usize offset = (row * mudGeom->getNumXCells() + column) * 3 + compIdx;
          equal &= (*mudRgb)[offset] == (*countRgb)[offset];
        }
      }
    }
    REQUIRE(equal);
  };
  // The source-unit legend differs, but all section and color-bar pixels must match.
  for(usize sectionIndex = 0; sectionIndex < 6; ++sectionIndex)
  {
    const auto origin = ebsdlib::GetODFSectionPanelOrigin(layout, sectionIndex);
    requireEqualRegion(static_cast<usize>(origin[0]), static_cast<usize>(origin[1]), static_cast<usize>(layout.panelWidth), static_cast<usize>(layout.panelHeight));
  }
  requireEqualRegion(static_cast<usize>(layout.columns * layout.panelSlotWidth + 2 * layout.margin), static_cast<usize>(layout.titleHeight + layout.fontPtSize + layout.margin),
                     static_cast<usize>(2 * layout.margin), static_cast<usize>(layout.panelHeight));
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Output modes preserve PNG and RGB pixels", "[MTRSim][WriteODFSectionsFilter][Execution]")
{
  for(const auto& mode : std::vector<std::pair<bool, bool>>{{false, true}, {true, false}, {true, true}})
  {
    DYNAMIC_SECTION("PNG " << mode.first << ", geometry " << mode.second)
    {
      TemporaryOutput output;
      auto data = CreateValidDataStructure();
      FillOdf(data, 1.0);
      auto args = ExecutionArguments();
      args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, mode.first);
      args.insertOrAssign(WriteODFSectionsFilter::k_SaveAsImageGeometry_Key, mode.second);
      args.insertOrAssign(WriteODFSectionsFilter::k_OutputPngPath_Key, output.png);
      const auto result = WriteODFSectionsFilter{}.execute(data, args);
      SIMPLNX_RESULT_REQUIRE_VALID(result.result);
      REQUIRE(std::filesystem::exists(output.png) == mode.first);
      const auto* geom = data.getDataAs<ImageGeom>(k_OutputGeometryPath);
      const auto* rgb = data.getDataAs<UInt8Array>(k_RgbPath);
      REQUIRE((geom != nullptr) == mode.second);
      REQUIRE((rgb != nullptr) == mode.second);
      if(mode.second)
      {
        REQUIRE(rgb->getNumberOfComponents() == 3);
        REQUIRE(std::any_of(rgb->begin(), rgb->end(), [](uint8 value) { return value != 0; }));
      }
      if(mode.first)
      {
        auto decoded = ReadPng(output.png);
        const auto* decodedGeom = decoded.getDataAs<ImageGeom>(DataPath({"Decoded"}));
        const auto* rgba = decoded.getDataAs<UInt8Array>(DataPath({"Decoded", "Cells", "Pixels"}));
        REQUIRE(decodedGeom != nullptr);
        REQUIRE(rgba != nullptr);
        REQUIRE(rgba->getNumberOfComponents() == 4);
        if(mode.second)
        {
          REQUIRE(decodedGeom->getDimensions() == geom->getDimensions());
          REQUIRE(rgba->getNumberOfTuples() == rgb->getNumberOfTuples());
          bool pixelsMatch = true;
          for(usize tupleIdx = 0; tupleIdx < rgb->getNumberOfTuples(); ++tupleIdx)
          {
            for(usize compIdx = 0; compIdx < 3; ++compIdx)
            {
              pixelsMatch &= (*rgb)[tupleIdx * 3 + compIdx] == (*rgba)[tupleIdx * 4 + compIdx];
            }
            pixelsMatch &= (*rgba)[tupleIdx * 4 + 3] == 255;
          }
          REQUIRE(pixelsMatch);
        }
      }
      REQUIRE(std::distance(std::filesystem::directory_iterator(output.directory), std::filesystem::directory_iterator{}) == (mode.first ? 1 : 0));
    }
  }
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Float32 and Float64 render identically", "[MTRSim][WriteODFSectionsFilter][Execution]")
{
  auto data32 = CreateValidDataStructure();
  auto data64 = CreateValidDataStructure();
  auto* cells = data32.getDataAs<AttributeMatrix>(k_OdfArrayPath.getParent());
  REQUIRE(cells != nullptr);
  auto* odf32 = Float32Array::CreateWithStore<DataStore<float32>>(data32, "Float ODF", {k_NumPhi1, k_NumPhi, k_NumPhi2}, {1}, cells->getId());
  auto* odf64 = data64.getDataAs<Float64Array>(k_OdfArrayPath);
  REQUIRE(odf32 != nullptr);
  REQUIRE(odf64 != nullptr);
  for(usize index = 0; index < odf32->getSize(); ++index)
  {
    (*odf32)[index] = 0.5F + static_cast<float>((index / k_NumPhi2) % 9) * 0.125F;
    (*odf64)[index] = (*odf32)[index];
  }
  auto args = ExecutionArguments();
  const auto result64 = WriteODFSectionsFilter{}.execute(data64, args);
  SIMPLNX_RESULT_REQUIRE_VALID(result64.result);
  args.insertOrAssign(WriteODFSectionsFilter::k_ODFArrayPath_Key, k_OdfArrayPath.replaceName("Float ODF"));
  const auto result32 = WriteODFSectionsFilter{}.execute(data32, args);
  SIMPLNX_RESULT_REQUIRE_VALID(result32.result);
  const auto* rgb32 = data32.getDataAs<UInt8Array>(k_RgbPath);
  const auto* rgb64 = data64.getDataAs<UInt8Array>(k_RgbPath);
  REQUIRE(rgb32 != nullptr);
  REQUIRE(rgb64 != nullptr);
  REQUIRE(std::equal(rgb32->begin(), rgb32->end(), rgb64->begin(), rgb64->end()));
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Invalid source values preserve an existing destination", "[MTRSim][WriteODFSectionsFilter][Execution]")
{
  const std::vector<std::pair<double, std::string>> cases = {{-2.5, "-2.5"}, {std::numeric_limits<double>::quiet_NaN(), "nan"}, {std::numeric_limits<double>::infinity(), "inf"}};
  for(const auto& [value, diagnostic] : cases)
  {
    DYNAMIC_SECTION(diagnostic)
    {
      TemporaryOutput output;
      const std::string original = "Existing destination must survive a failed render.\n";
      {
        std::ofstream file(output.png, std::ios::binary);
        file << original;
        REQUIRE(file.good());
      }
      auto data = CreateValidDataStructure();
      FillOdf(data, 1.0);
      auto* odf = data.getDataAs<Float64Array>(k_OdfArrayPath);
      REQUIRE(odf != nullptr);
      // This index is outside the displayed crop and still requires validation.
      (*odf)[odf->getSize() - 1] = value;
      auto args = ExecutionArguments();
      args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, true);
      args.insertOrAssign(WriteODFSectionsFilter::k_OutputPngPath_Key, output.png);
      const auto result = WriteODFSectionsFilter{}.execute(data, args);
      SIMPLNX_RESULT_REQUIRE_INVALID(result.result);
      const auto& message = result.result.errors().front().message;
      REQUIRE(message.find(k_OdfArrayPath.toString()) != std::string::npos);
      REQUIRE(message.find("phase 1") != std::string::npos);
      REQUIRE(message.find("index 186623") != std::string::npos);
      REQUIRE(message.find(diagnostic) != std::string::npos);
      std::ifstream file(output.png, std::ios::binary);
      REQUIRE(std::string(std::istreambuf_iterator<char>(file), {}) == original);
      REQUIRE(std::distance(std::filesystem::directory_iterator(output.directory), std::filesystem::directory_iterator{}) == 1);
    }
  }
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Zero data requires manual scaling and uses the first preset color", "[MTRSim][WriteODFSectionsFilter][Execution]")
{
  auto data = CreateValidDataStructure();
  FillOdf(data, 0.0);
  auto args = ExecutionArguments();
  SECTION("Automatic")
  {
    const auto result = WriteODFSectionsFilter{}.execute(data, args);
    SIMPLNX_RESULT_REQUIRE_INVALID(result.result);
    REQUIRE(result.result.errors().front().message.find("maximum is 0") != std::string::npos);
  }
  SECTION("Manual")
  {
    const auto preset = GENERATE(std::string{"Rainbow"}, std::string{"Grayscale"});
    args.insertOrAssign(WriteODFSectionsFilter::k_ColorTablePreset_Key, preset);
    args.insertOrAssign(WriteODFSectionsFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
    args.insertOrAssign(WriteODFSectionsFilter::k_ManualMinimum_Key, float64{0.25});
    args.insertOrAssign(WriteODFSectionsFilter::k_ManualMaximum_Key, float64{2.5});
    const auto result = WriteODFSectionsFilter{}.execute(data, args);
    SIMPLNX_RESULT_REQUIRE_VALID(result.result);
    const std::array<uint8, 3> firstColor = preset == "Rainbow" ? std::array<uint8, 3>{0, 0, 255} : std::array<uint8, 3>{0, 0, 0};
    ebsdlib::ODFSectionConfiguration config;
    config.grid.laueOpsIndex = ebsdlib::CrystalStructure::Hexagonal_High;
    config.sectionWidth = 128;
    const auto layout = ebsdlib::ComputeODFSectionLayout(config);
    const auto* rgb = data.getDataAs<UInt8Array>(k_RgbPath);
    REQUIRE(rgb != nullptr);
    for(usize sectionIndex = 0; sectionIndex < 6; ++sectionIndex)
    {
      const auto origin = ebsdlib::GetODFSectionPanelOrigin(layout, sectionIndex);
      const usize pixelIndex = (static_cast<usize>(origin[1]) + static_cast<usize>(layout.panelHeight / 2)) * layout.pageWidth + static_cast<usize>(origin[0]) + 64;
      for(usize compIdx = 0; compIdx < 3; ++compIdx)
      {
        REQUIRE((*rgb)[pixelIndex * 3 + compIdx] == firstColor[compIdx]);
      }
    }
  }
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Cancellation at each stage creates no destination", "[MTRSim][WriteODFSectionsFilter][Execution]")
{
  const std::vector<std::string> stages = {"Converting ODF values", "Composing ODF sections", "Copying ODF section RGB pixels", "Writing ODF section PNG"};
  for(const auto& stage : stages)
  {
    DYNAMIC_SECTION(stage)
    {
      TemporaryOutput output;
      auto data = CreateValidDataStructure();
      FillOdf(data, 1.0);
      auto args = ExecutionArguments();
      args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, true);
      args.insertOrAssign(WriteODFSectionsFilter::k_OutputPngPath_Key, output.png);
      std::atomic_bool cancel = false;
      bool reachedStage = false;
      const IFilter::MessageHandler handler{[&](const IFilter::Message& message) {
        if(message.message == "Converting ODF values")
        {
          auto* rgb = data.getDataAs<UInt8Array>(k_RgbPath);
          REQUIRE(rgb != nullptr);
          std::fill(rgb->begin(), rgb->end(), uint8{173});
        }
        if(message.message == stage)
        {
          reachedStage = true;
          cancel = true;
        }
      }};
      const auto result = WriteODFSectionsFilter{}.execute(data, args, nullptr, handler, cancel);
      SIMPLNX_RESULT_REQUIRE_INVALID(result.result);
      REQUIRE(result.result.errors().front().code == -1);
      REQUIRE(reachedStage);
      REQUIRE_FALSE(std::filesystem::exists(output.png));
      REQUIRE(std::filesystem::is_empty(output.directory));
      if(stage != "Writing ODF section PNG")
      {
        const auto* rgb = data.getDataAs<UInt8Array>(k_RgbPath);
        REQUIRE(rgb != nullptr);
        REQUIRE(std::all_of(rgb->begin(), rgb->end(), [](uint8 value) { return value == 173; }));
      }
    }
  }
  SECTION("Already cancelled algorithm does not access source data")
  {
    DataStructure data;
    WriteODFSectionsInputValues values;
    const std::atomic_bool cancel = true;
    const IFilter::MessageHandler handler;
    const auto result = WriteODFSections(data, handler, cancel, &values)();
    SIMPLNX_RESULT_REQUIRE_VALID(result);
    REQUIRE(data.getSize() == 0);
  }
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Compute ODF produces six HCP sections and a readable PNG", "[MTRSim][WriteODFSectionsFilter][Execution]")
{
  const float32 binSize = GENERATE(5.0F, static_cast<float32>(180.0 / 7.0));
  TemporaryOutput output;
  auto data = CreateValidDataStructure();
  auto* cells = AttributeMatrix::Create(data, "Input Cells", {8});
  REQUIRE(cells != nullptr);
  auto* eulers = Float32Array::CreateWithStore<DataStore<float32>>(data, "Euler Angles", {8}, {3}, cells->getId());
  auto* phases = Int32Array::CreateWithStore<DataStore<int32>>(data, "Phases", {8}, {1}, cells->getId());
  REQUIRE(eulers != nullptr);
  REQUIRE(phases != nullptr);
  for(usize tupleIdx = 0; tupleIdx < 8; ++tupleIdx)
  {
    (*eulers)[tupleIdx * 3] = static_cast<float32>((15.0 + tupleIdx * 10.0) * std::numbers::pi / 180.0);
    (*eulers)[tupleIdx * 3 + 1] = static_cast<float32>((20.0 + tupleIdx * 5.0) * std::numbers::pi / 180.0);
    (*eulers)[tupleIdx * 3 + 2] = static_cast<float32>((5.0 + tupleIdx * 5.0) * std::numbers::pi / 180.0);
    (*phases)[tupleIdx] = 1;
  }
  ComputeODFFilter computeFilter;
  auto computeArgs = computeFilter.getDefaultArguments();
  computeArgs.insertOrAssign(ComputeODFFilter::k_EulerAngles_Key, DataPath({"Input Cells", "Euler Angles"}));
  computeArgs.insertOrAssign(ComputeODFFilter::k_Phases_Key, DataPath({"Input Cells", "Phases"}));
  computeArgs.insertOrAssign(ComputeODFFilter::k_CrystalStructures_Key, k_CrystalStructuresPath);
  computeArgs.insertOrAssign(ComputeODFFilter::k_UseMask_Key, false);
  computeArgs.insertOrAssign(ComputeODFFilter::k_ApplySmoothing_Key, true);
  computeArgs.insertOrAssign(ComputeODFFilter::k_BinSizeDeg_Key, binSize);
  computeArgs.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, ChoicesParameter::ValueType{1});
  computeArgs.insertOrAssign(ComputeODFFilter::k_OutputImageGeometry_Key, DataPath({"Computed ODF"}));
  computeArgs.insertOrAssign(ComputeODFFilter::k_CellAttrMatName_Key, std::string{"Cells"});
  computeArgs.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::string{"MUD"});
  const auto computePreflight = computeFilter.preflight(data, computeArgs);
  SIMPLNX_RESULT_REQUIRE_VALID(computePreflight.outputActions);
  const auto computed = computeFilter.execute(data, computeArgs);
  SIMPLNX_RESULT_REQUIRE_VALID(computed.result);
  auto args = ExecutionArguments();
  args.insertOrAssign(WriteODFSectionsFilter::k_InputImageGeometryPath_Key, DataPath({"Computed ODF"}));
  args.insertOrAssign(WriteODFSectionsFilter::k_ODFArrayPath_Key, DataPath({"Computed ODF", "Cells", "MUD"}));
  args.insertOrAssign(WriteODFSectionsFilter::k_WritePng_Key, true);
  args.insertOrAssign(WriteODFSectionsFilter::k_OutputPngPath_Key, output.png);
  WriteODFSectionsFilter filter;
  const auto preflight = filter.preflight(data, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.outputActions);
  REQUIRE(GetPreflightValue(preflight, "Section Angles") == "0, 10, 20, 30, 40, 50 degrees");
  double maximum = std::numeric_limits<double>::quiet_NaN();
  const IFilter::MessageHandler handler{[&](const IFilter::Message& message) {
    const std::string marker = "maximum MUD = ";
    const auto position = message.message.find(marker);
    if(position != std::string::npos)
    {
      maximum = std::stod(message.message.substr(position + marker.size()));
    }
  }};
  const auto result = filter.execute(data, args, nullptr, handler);
  SIMPLNX_RESULT_REQUIRE_VALID(result.result);
  REQUIRE(std::isfinite(maximum));
  REQUIRE(maximum > 0.0);
  const auto* rgb = data.getDataAs<UInt8Array>(k_RgbPath);
  REQUIRE(rgb != nullptr);
  REQUIRE(rgb->getNumberOfTuples() > 0);
  auto decoded = ReadPng(output.png);
  REQUIRE(decoded.getDataAs<ImageGeom>(DataPath({"Decoded"})) != nullptr);
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Deferred metadata preflight uses a conservative layout", "[MTRSim][WriteODFSectionsFilter][Deferred]")
{
  auto data = CreateDeferredMetadataDataStructure();
  auto* crystals = data.getDataAs<UInt32Array>(k_CrystalStructuresPath);
  REQUIRE(crystals != nullptr);
  REQUIRE(crystals->getIDataStoreRef().getStoreType() == IDataStore::StoreType::Empty);
  SECTION("Empty in-memory store")
  {
  }
  SECTION("Empty out-of-core store")
  {
    crystals->setDataStore(std::make_shared<EmptyDataStore<uint32>>(ShapeType{2}, ShapeType{1}, "deferred-test"));
    REQUIRE(crystals->getIDataStoreRef().getStoreType() == IDataStore::StoreType::EmptyOutOfCore);
  }
  IFilter::PreflightResult preflight;
  REQUIRE_NOTHROW(preflight = WriteODFSectionsFilter{}.preflight(data, ExecutionArguments()));
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.outputActions);
  for(const auto& label : {"Selected Phase Metadata", "Laue Group", "Section Angles"})
  {
    REQUIRE(GetPreflightValue(preflight, label).find("Deferred") != std::string::npos);
  }
  REQUIRE(GetPreflightValue(preflight, "Composite Dimensions").find("conservative maximum") != std::string::npos);
  REQUIRE(GetPreflightValue(preflight, "Input Units") == "MUD");
  REQUIRE(GetPreflightValue(preflight, "Section Layout") == "3 columns x 2 rows");
  REQUIRE(preflight.outputActions.value().actions.size() == 2);
  const auto* geomAction = dynamic_cast<const CreateImageGeometryAction*>(preflight.outputActions.value().actions[0].get());
  const auto* rgbAction = dynamic_cast<const CreateArrayAction*>(preflight.outputActions.value().actions[1].get());
  REQUIRE(geomAction != nullptr);
  REQUIRE(rgbAction != nullptr);
  // Width 128, Triclinic panel height 64, and two panel rows give a 628 by 270 page.
  REQUIRE(geomAction->dims() == CreateImageGeometryAction::DimensionType{628, 270, 1});
  REQUIRE(rgbAction->dims() == ShapeType{1, 270, 628});
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Deferred metadata execution shrinks the conservative geometry", "[MTRSim][WriteODFSectionsFilter][Deferred]")
{
  TemporaryOutput output;
  auto data = CreateDeferredMetadataDataStructure();
  FillOdf(data, 1.0);
  IFilter::PreflightResult preflight;
  REQUIRE_NOTHROW(preflight = WriteODFSectionsFilter{}.preflight(data, ExecutionArguments()));
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.outputActions);
  const auto applied = preflight.outputActions.value().applyRegular(data, IDataAction::Mode::Execute);
  SIMPLNX_RESULT_REQUIRE_VALID(applied);
  const auto* geometry = data.getDataAs<ImageGeom>(k_OutputGeometryPath);
  REQUIRE(geometry != nullptr);
  REQUIRE(geometry->getDimensions() == SizeVec3{628, 270, 1});
  auto* crystals = data.getDataAs<UInt32Array>(k_CrystalStructuresPath);
  REQUIRE(crystals != nullptr);
  crystals->setDataStore(std::make_shared<DataStore<uint32>>(ShapeType{2}, ShapeType{1}, ebsdlib::CrystalStructure::UnknownCrystalStructure));
  (*crystals)[1] = ebsdlib::CrystalStructure::Hexagonal_High;
  auto* names = data.getDataAs<StringArray>(k_MaterialNamesPath);
  REQUIRE(names != nullptr);
  (*names)[1] = "Titanium (Alpha)";
  auto values = ExecutionInputValues();
  values.writePng = true;
  values.outputPngPath = output.png;
  const std::atomic_bool cancel = false;
  const IFilter::MessageHandler handler;
  const auto result = WriteODFSections(data, handler, cancel, &values)();
  SIMPLNX_RESULT_REQUIRE_VALID(result);
  REQUIRE(geometry->getDimensions() == SizeVec3{628, 206, 1});
  REQUIRE(geometry->getCellData()->getShape() == ShapeType{1, 206, 628});
  const auto* rgb = data.getDataAs<UInt8Array>(k_RgbPath);
  REQUIRE(rgb != nullptr);
  REQUIRE(rgb->getTupleShape() == ShapeType{1, 206, 628});
  auto decoded = ReadPng(output.png);
  const auto* decodedGeometry = decoded.getDataAs<ImageGeom>(DataPath({"Decoded"}));
  const auto* rgba = decoded.getDataAs<UInt8Array>(DataPath({"Decoded", "Cells", "Pixels"}));
  REQUIRE(decodedGeometry != nullptr);
  REQUIRE(rgba != nullptr);
  REQUIRE(decodedGeometry->getDimensions() == geometry->getDimensions());
  bool matches = true;
  for(usize tupleIdx = 0; tupleIdx < rgb->getNumberOfTuples(); ++tupleIdx)
  {
    for(usize compIdx = 0; compIdx < 3; ++compIdx)
    {
      matches &= (*rgb)[tupleIdx * 3 + compIdx] == (*rgba)[tupleIdx * 4 + compIdx];
    }
  }
  REQUIRE(matches);
}

TEST_CASE("MTRSim::WriteODFSectionsFilter: Deferred metadata execution rejects missing or invalid phase values", "[MTRSim][WriteODFSectionsFilter][Deferred]")
{
  auto data = CreateDeferredMetadataDataStructure();
  FillOdf(data, 1.0);
  auto values = ExecutionInputValues();
  values.saveAsImageGeometry = false;
  auto* crystals = data.getDataAs<UInt32Array>(k_CrystalStructuresPath);
  REQUIRE(crystals != nullptr);
  std::string expectedFragment;
  SECTION("Still empty")
  {
    expectedFragment = "not available";
  }
  SECTION("Still empty out-of-core")
  {
    crystals->setDataStore(std::make_shared<EmptyDataStore<uint32>>(ShapeType{2}, ShapeType{1}, "deferred-test"));
    expectedFragment = "not available";
  }
  SECTION("Unknown crystal structure after materialization")
  {
    crystals->setDataStore(std::make_shared<DataStore<uint32>>(ShapeType{2}, ShapeType{1}, uint32{999}));
    expectedFragment = "999";
  }
  SECTION("Phase zero after materialization")
  {
    crystals->setDataStore(std::make_shared<DataStore<uint32>>(ShapeType{2}, ShapeType{1}, ebsdlib::CrystalStructure::Hexagonal_High));
    values.phaseIndex = 0;
    expectedFragment = "Phase Index (0)";
  }
  SECTION("Phase outside materialized tuple range")
  {
    crystals->setDataStore(std::make_shared<DataStore<uint32>>(ShapeType{2}, ShapeType{1}, ebsdlib::CrystalStructure::Hexagonal_High));
    values.phaseIndex = 2;
    expectedFragment = "Phase Index (2)";
  }
  SECTION("Materialized ensemble tuple counts differ")
  {
    crystals->setDataStore(std::make_shared<DataStore<uint32>>(ShapeType{3}, ShapeType{1}, ebsdlib::CrystalStructure::Hexagonal_High));
    expectedFragment = "contains 3 tuples";
  }
  const std::atomic_bool cancel = false;
  const IFilter::MessageHandler handler;
  Result<> result;
  REQUIRE_NOTHROW(result = WriteODFSections(data, handler, cancel, &values)());
  SIMPLNX_RESULT_REQUIRE_INVALID(result);
  REQUIRE(result.errors().front().code <= -12300);
  REQUIRE(result.errors().front().code > -12400);
  REQUIRE(result.errors().front().message.find(k_CrystalStructuresPath.toString()) != std::string::npos);
  REQUIRE(result.errors().front().message.find(expectedFragment) != std::string::npos);
}

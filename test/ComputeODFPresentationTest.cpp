#include <catch2/catch.hpp>

#include "MTRSim/Filters/Algorithms/GenerateODFSections.hpp"
#include "MTRSim/Filters/ComputeODFFilter.hpp"

#include "NXImageProcessing/Filters/ReadImageFilter.hpp"

#include "EbsdLib/Core/EbsdLibConstants.h"
#include "EbsdLib/Utilities/ODFSectionChrome.h"
#include "EbsdLib/Utilities/ODFSectionCompositor.h"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataStore.hpp"
#include "simplnx/DataStructure/EmptyStringStore.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/StringArray.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/CreateColorMapParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"

#include <algorithm>
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
const DataPath k_CrystalStructuresPath({"Ensemble Data", "Crystal Structures"});
const DataPath k_MaterialNamesPath({"Ensemble Data", "Material Names"});
const DataPath k_OdfGeometryPath({"ODF"});
const DataPath k_OdfArrayPath({"ODF", "Cell Data", "ODF Values"});
const DataPath k_OutputGeometryPath({"ODF Sections"});
const DataPath k_RgbPath({"ODF Sections", "Cell Data", "ODF Sections RGB"});

DataStructure CreateValidDataStructure()
{
  DataStructure data;
  auto* cellsPtr = AttributeMatrix::Create(data, "Input Cells", {8});
  auto* eulersPtr = Float32Array::Create(data, "Euler Angles", DataStoreUtilities::CreateDataStore<float32>(data, DataPath({"Input Cells", "Euler Angles"}), {8}, {3}), cellsPtr->getId());
  auto* phasesPtr = Int32Array::Create(data, "Phases", DataStoreUtilities::CreateDataStore<int32>(data, DataPath({"Input Cells", "Phases"}), {8}, {1}), cellsPtr->getId());
  REQUIRE(eulersPtr != nullptr);
  REQUIRE(phasesPtr != nullptr);
  for(usize tupleIdx = 0; tupleIdx < 8; ++tupleIdx)
  {
    (*eulersPtr)[tupleIdx * 3] = static_cast<float32>((15.0 + tupleIdx * 10.0) * std::numbers::pi / 180.0);
    (*eulersPtr)[tupleIdx * 3 + 1] = static_cast<float32>((20.0 + tupleIdx * 5.0) * std::numbers::pi / 180.0);
    (*eulersPtr)[tupleIdx * 3 + 2] = static_cast<float32>((5.0 + tupleIdx * 5.0) * std::numbers::pi / 180.0);
    (*phasesPtr)[tupleIdx] = 1;
  }
  auto* ensemblePtr = AttributeMatrix::Create(data, "Ensemble Data", {2});
  auto* crystalsPtr = UInt32Array::Create(data, "Crystal Structures", DataStoreUtilities::CreateDataStore<uint32>(data, k_CrystalStructuresPath, {2}, {1}), ensemblePtr->getId());
  REQUIRE(crystalsPtr != nullptr);
  (*crystalsPtr)[0] = ebsdlib::CrystalStructure::UnknownCrystalStructure;
  (*crystalsPtr)[1] = ebsdlib::CrystalStructure::Hexagonal_High;
  REQUIRE(StringArray::CreateWithValues(data, "Material Names", {2}, {"Invalid", "Titanium (Alpha)"}, ensemblePtr->getId()) != nullptr);
  return data;
}

Arguments ExecutionArguments()
{
  auto args = ComputeODFFilter{}.getDefaultArguments();
  args.insertOrAssign(ComputeODFFilter::k_EulerAngles_Key, DataPath({"Input Cells", "Euler Angles"}));
  args.insertOrAssign(ComputeODFFilter::k_Phases_Key, DataPath({"Input Cells", "Phases"}));
  args.insertOrAssign(ComputeODFFilter::k_CrystalStructures_Key, k_CrystalStructuresPath);
  args.insertOrAssign(ComputeODFFilter::k_MaterialNames_Key, k_MaterialNamesPath);
  args.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::string{"ODF Values"});
  args.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, ChoicesParameter::ValueType{1});
  args.insertOrAssign("section_width", int32{128});
  args.insertOrAssign("write_png", false);
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
  int32 error_code = 0;
  std::string diagnostic_fragment;
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
  const auto* ensemblePtr = data.getDataAs<AttributeMatrix>(k_CrystalStructuresPath.getParent());
  REQUIRE(ensemblePtr != nullptr);
  REQUIRE(StringArray::CreateWithValues(data, k_MaterialNamesPath.getTargetName(), {2}, {"", ""}, ensemblePtr->getId()) != nullptr);
  return data;
}

GenerateODFSectionsInputValues ExecutionInputValues()
{
  GenerateODFSectionsInputValues values;
  values.inputImageGeometryPath = k_OdfGeometryPath;
  values.odfArrayPath = k_OdfArrayPath;
  values.crystalStructure = ebsdlib::CrystalStructure::Hexagonal_High;
  values.materialName = "Titanium (Alpha)";
  values.title = "ODF Sections";
  values.sectionWidth = 128;
  return values;
}

DataStructure ComputeNumericData()
{
  auto data = CreateValidDataStructure();
  auto args = ExecutionArguments();
  args.insertOrAssign(ComputeODFFilter::k_SaveAsImageGeometry_Key, false);
  const auto result = ComputeODFFilter{}.execute(data, args);
  SIMPLNX_RESULT_REQUIRE_VALID(result.result);
  return data;
}

Result<> RenderExisting(DataStructure& data, const GenerateODFSectionsInputValues& values)
{
  if(values.saveAsImageGeometry)
  {
    ebsdlib::ODFSectionConfiguration config;
    config.grid.laueOpsIndex = values.crystalStructure;
    config.sectionCount = static_cast<usize>(values.numberOfSections);
    config.sectionsPerRow = static_cast<usize>(values.sectionsPerRow);
    config.sectionWidth = values.sectionWidth;
    const auto layout = ebsdlib::ComputeODFSectionLayout(config);
    const CreateImageGeometryAction geometryAction(values.outputImageGeometryPath, {static_cast<usize>(layout.pageWidth), static_cast<usize>(layout.pageHeight), 1}, {0, 0, 0}, {1, 1, 1}, "Cell Data");
    const auto geometryCreated = geometryAction.apply(data, IDataAction::Mode::Execute);
    SIMPLNX_RESULT_REQUIRE_VALID(geometryCreated);
    const CreateArrayAction arrayAction(DataType::uint8, {1, static_cast<usize>(layout.pageHeight), static_cast<usize>(layout.pageWidth)}, {3}, k_RgbPath);
    const auto arrayCreated = arrayAction.apply(data, IDataAction::Mode::Execute);
    SIMPLNX_RESULT_REQUIRE_VALID(arrayCreated);
  }
  const IFilter::MessageHandler handler;
  const std::atomic_bool cancel = false;
  return GenerateODFSections(data, handler, cancel, &values)();
}
} // namespace

TEST_CASE("MTRSim::ComputeODFFilter: Output modes preserve PNG and RGB pixels", "[MTRSim][ComputeODFFilter][Execution]")
{
  for(const auto& mode : std::vector<std::pair<bool, bool>>{{false, false}, {false, true}, {true, false}, {true, true}})
  {
    DYNAMIC_SECTION("PNG " << mode.first << ", geometry " << mode.second)
    {
      TemporaryOutput output;
      auto data = CreateValidDataStructure();
      auto args = ExecutionArguments();
      args.insertOrAssign("write_png", mode.first);
      args.insertOrAssign("save_as_image_geometry", mode.second);
      args.insertOrAssign("output_png_file", output.png);
      const auto result = ComputeODFFilter{}.execute(data, args);
      SIMPLNX_RESULT_REQUIRE_VALID(result.result);
      REQUIRE(std::filesystem::exists(output.png) == mode.first);
      const auto* geomPtr = data.getDataAs<ImageGeom>(k_OutputGeometryPath);
      const auto* rgbPtr = data.getDataAs<UInt8Array>(k_RgbPath);
      REQUIRE((geomPtr != nullptr) == mode.second);
      REQUIRE((rgbPtr != nullptr) == mode.second);
      if(mode.second)
      {
        REQUIRE(rgbPtr->getNumberOfComponents() == 3);
        REQUIRE(std::any_of(rgbPtr->begin(), rgbPtr->end(), [](uint8 value) { return value != 0; }));
      }
      if(mode.first)
      {
        auto decoded = ReadPng(output.png);
        const auto* decodedGeomPtr = decoded.getDataAs<ImageGeom>(DataPath({"Decoded"}));
        const auto* rgbaPtr = decoded.getDataAs<UInt8Array>(DataPath({"Decoded", "Cells", "Pixels"}));
        REQUIRE(decodedGeomPtr != nullptr);
        REQUIRE(rgbaPtr != nullptr);
        REQUIRE(rgbaPtr->getNumberOfComponents() == 4);
        if(mode.second)
        {
          REQUIRE(decodedGeomPtr->getDimensions() == geomPtr->getDimensions());
          REQUIRE(rgbaPtr->getNumberOfTuples() == rgbPtr->getNumberOfTuples());
          bool pixelsMatch = true;
          for(usize tupleIdx = 0; tupleIdx < rgbPtr->getNumberOfTuples(); ++tupleIdx)
          {
            for(usize compIdx = 0; compIdx < 3; ++compIdx)
            {
              pixelsMatch &= (*rgbPtr)[tupleIdx * 3 + compIdx] == (*rgbaPtr)[tupleIdx * 4 + compIdx];
            }
            pixelsMatch &= (*rgbaPtr)[tupleIdx * 4 + 3] == 255;
          }
          REQUIRE(pixelsMatch);
        }
      }
      REQUIRE(std::distance(std::filesystem::directory_iterator(output.directory), std::filesystem::directory_iterator{}) == (mode.first ? 1 : 0));
    }
  }
}

TEST_CASE("MTRSim::ComputeODFFilter: Presentation preflight supplies exact phase metadata", "[MTRSim][ComputeODFFilter][Presentation]")
{
  auto data = CreateValidDataStructure();
  const auto result = ComputeODFFilter{}.preflight(data, ExecutionArguments());
  SIMPLNX_RESULT_REQUIRE_VALID(result.outputActions);
  REQUIRE(GetPreflightValue(result, "Section Angles") == "0, 10, 20, 30, 40, 50 degrees");
  REQUIRE(GetPreflightValue(result, "Material") == "Titanium (Alpha)");
  REQUIRE(GetPreflightValue(result, "Laue Class") == "Hexagonal 6/mmm (D6h)");
  REQUIRE(GetPreflightValue(result, "Rendered Units") == "MUD");
}
TEST_CASE("MTRSim::ComputeODFFilter: Units produce identical MUD and Count-Density section pixels", "[MTRSim][ComputeODFFilter][Units]")
{
  auto mudData = CreateValidDataStructure();
  auto countData = CreateValidDataStructure();
  auto args = ExecutionArguments();
  args.insertOrAssign(ComputeODFFilter::k_WritePng_Key, false);
  args.insertOrAssign(ComputeODFFilter::k_SectionWidth_Key, int32{128});
  args.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, ChoicesParameter::ValueType{1});
  ComputeODFFilter filter;
  const auto mudResult = filter.execute(mudData, args);
  SIMPLNX_RESULT_REQUIRE_VALID(mudResult.result);
  args.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, ChoicesParameter::ValueType{0});
  const auto countResult = filter.execute(countData, args);
  SIMPLNX_RESULT_REQUIRE_VALID(countResult.result);
  const auto* mudGeomPtr = mudData.getDataAs<ImageGeom>(k_OutputGeometryPath);
  const auto* countGeomPtr = countData.getDataAs<ImageGeom>(k_OutputGeometryPath);
  REQUIRE(mudGeomPtr != nullptr);
  REQUIRE(countGeomPtr != nullptr);
  REQUIRE(mudGeomPtr->getDimensions() == countGeomPtr->getDimensions());
  const auto rgbPath = k_OutputGeometryPath.createChildPath("Cell Data").createChildPath("ODF Sections RGB");
  const auto* mudRgbPtr = mudData.getDataAs<UInt8Array>(rgbPath);
  const auto* countRgbPtr = countData.getDataAs<UInt8Array>(rgbPath);
  REQUIRE(mudRgbPtr != nullptr);
  REQUIRE(countRgbPtr != nullptr);
  REQUIRE(ebsdlib::GetODFColorBarTitle(ebsdlib::ODFValueUnits::CountDensity) == "MUD (from Count-Density)");
  REQUIRE_FALSE(std::equal(mudRgbPtr->begin(), mudRgbPtr->end(), countRgbPtr->begin(), countRgbPtr->end()));
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
          const usize offset = (row * mudGeomPtr->getNumXCells() + column) * 3 + compIdx;
          equal &= (*mudRgbPtr)[offset] == (*countRgbPtr)[offset];
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

TEST_CASE("MTRSim::ComputeODFFilter: Float32 fractional bin spacing reaches presentation", "[MTRSim][ComputeODFFilter][Presentation]")
{
  const usize phiCount = GENERATE(usize{7}, usize{53});
  auto data = CreateValidDataStructure();
  auto args = ExecutionArguments();
  const float32 spacing = static_cast<float32>(180.0 / static_cast<double>(phiCount));
  args.insertOrAssign(ComputeODFFilter::k_BinSizeDeg_Key, spacing);
  const auto result = ComputeODFFilter{}.execute(data, args);
  SIMPLNX_RESULT_REQUIRE_VALID(result.result);
  REQUIRE(data.getDataAs<UInt8Array>(k_RgbPath) != nullptr);
  auto invalidArgs = ExecutionArguments();
  invalidArgs.insertOrAssign(ComputeODFFilter::k_BinSizeDeg_Key, spacing * 1.000001F);
  auto invalidData = CreateValidDataStructure();
  const auto invalid = ComputeODFFilter{}.preflight(invalidData, invalidArgs);
  SIMPLNX_RESULT_REQUIRE_INVALID(invalid.outputActions);
  REQUIRE(invalid.outputActions.errors().front().code == -12200);
}
TEST_CASE("MTRSim::ComputeODFFilter: Cancellation at each stage creates no destination", "[MTRSim][ComputeODFFilter][Execution]")
{
  const std::vector<std::string> stages = {"Converting ODF values", "Composing ODF sections", "Copying ODF section RGB pixels", "Writing ODF section PNG"};
  for(const auto& stage : stages)
  {
    DYNAMIC_SECTION(stage)
    {
      TemporaryOutput output;
      auto data = CreateValidDataStructure();
      auto args = ExecutionArguments();
      args.insertOrAssign(ComputeODFFilter::k_WritePng_Key, true);
      args.insertOrAssign(ComputeODFFilter::k_OutputPngPath_Key, output.png);
      std::atomic_bool cancel = false;
      bool reachedStage = false;
      const IFilter::MessageHandler handler{[&](const IFilter::Message& message) {
        if(message.message == "Converting ODF values")
        {
          auto* rgbPtr = data.getDataAs<UInt8Array>(k_RgbPath);
          REQUIRE(rgbPtr != nullptr);
          std::fill(rgbPtr->begin(), rgbPtr->end(), uint8{173});
        }
        if(message.message == stage)
        {
          reachedStage = true;
          cancel = true;
        }
      }};
      const auto result = ComputeODFFilter{}.execute(data, args, nullptr, handler, cancel);
      SIMPLNX_RESULT_REQUIRE_INVALID(result.result);
      REQUIRE(result.result.errors().front().code == -1);
      REQUIRE(reachedStage);
      REQUIRE_FALSE(std::filesystem::exists(output.png));
      REQUIRE(std::filesystem::is_empty(output.directory));
      if(stage != "Writing ODF section PNG")
      {
        const auto* rgbPtr = data.getDataAs<UInt8Array>(k_RgbPath);
        REQUIRE(rgbPtr != nullptr);
        REQUIRE(std::all_of(rgbPtr->begin(), rgbPtr->end(), [](uint8 value) { return value == 173; }));
      }
    }
  }
  SECTION("Already cancelled algorithm does not access source data")
  {
    DataStructure data;
    GenerateODFSectionsInputValues values;
    const std::atomic_bool cancel = true;
    const IFilter::MessageHandler handler;
    const auto result = GenerateODFSections(data, handler, cancel, &values)();
    SIMPLNX_RESULT_REQUIRE_VALID(result);
    REQUIRE(data.getSize() == 0);
  }
}

TEST_CASE("MTRSim::ComputeODFFilter: Presentation preflight rejects invalid settings", "[MTRSim][ComputeODFFilter][Presentation]")
{
  const std::vector<InvalidCase> cases = {
      {.name = "section count",
       .error_code = -12309,
       .diagnostic_fragment = "1",
       .mutate = [](DataStructure&, Arguments& args) { args.insertOrAssign(ComputeODFFilter::k_NumSections_Key, int32{1}); }},
      {.name = "sections per row zero",
       .error_code = -12310,
       .diagnostic_fragment = "0",
       .mutate = [](DataStructure&, Arguments& args) { args.insertOrAssign(ComputeODFFilter::k_SectionsPerRow_Key, int32{0}); }},
      {.name = "sections per row too large",
       .error_code = -12310,
       .diagnostic_fragment = "7",
       .mutate = [](DataStructure&, Arguments& args) { args.insertOrAssign(ComputeODFFilter::k_SectionsPerRow_Key, int32{7}); }},
      {.name = "section width",
       .error_code = -12311,
       .diagnostic_fragment = "127",
       .mutate = [](DataStructure&, Arguments& args) { args.insertOrAssign(ComputeODFFilter::k_SectionWidth_Key, int32{127}); }},
      {.name = "invalid preset",
       .error_code = -12312,
       .diagnostic_fragment = "Not A Preset",
       .mutate = [](DataStructure&, Arguments& args) { args.insertOrAssign(ComputeODFFilter::k_ColorTablePreset_Key, CreateColorMapParameter::ValueType{"Not A Preset"}); }},
      {.name = "negative manual minimum",
       .error_code = -12313,
       .diagnostic_fragment = "-0.5",
       .mutate =
           [](DataStructure&, Arguments& args) {
             args.insertOrAssign(ComputeODFFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
             args.insertOrAssign(ComputeODFFilter::k_ManualMinimum_Key, float64{-0.5});
           }},
      {.name = "non-finite manual minimum",
       .error_code = -12313,
       .diagnostic_fragment = "nan",
       .mutate =
           [](DataStructure&, Arguments& args) {
             args.insertOrAssign(ComputeODFFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
             args.insertOrAssign(ComputeODFFilter::k_ManualMinimum_Key, std::numeric_limits<float64>::quiet_NaN());
           }},
      {.name = "equal manual limits",
       .error_code = -12313,
       .diagnostic_fragment = "[1, 1]",
       .mutate =
           [](DataStructure&, Arguments& args) {
             args.insertOrAssign(ComputeODFFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
             args.insertOrAssign(ComputeODFFilter::k_ManualMinimum_Key, float64{1.0});
             args.insertOrAssign(ComputeODFFilter::k_ManualMaximum_Key, float64{1.0});
           }},
      {.name = "reversed manual limits",
       .error_code = -12313,
       .diagnostic_fragment = "[2, 1]",
       .mutate =
           [](DataStructure&, Arguments& args) {
             args.insertOrAssign(ComputeODFFilter::k_ScaleMode_Key, ChoicesParameter::ValueType{1});
             args.insertOrAssign(ComputeODFFilter::k_ManualMinimum_Key, float64{2.0});
             args.insertOrAssign(ComputeODFFilter::k_ManualMaximum_Key, float64{1.0});
           }},
      {.name = "non-PNG extension",
       .error_code = -3003,
       .diagnostic_fragment = ".jpg",
       .mutate =
           [](DataStructure&, Arguments& args) {
             args.insertOrAssign(ComputeODFFilter::k_WritePng_Key, true);
             args.insertOrAssign(ComputeODFFilter::k_OutputPngPath_Key, FileSystemPathParameter::ValueType{"odf_sections.jpg"});
           }},
  };

  for(const auto& testCase : cases)
  {
    DYNAMIC_SECTION(testCase.name)
    {
      DataStructure dataStructure = CreateValidDataStructure();
      Arguments args = ExecutionArguments();
      testCase.mutate(dataStructure, args);

      ComputeODFFilter filter;
      const auto result = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_INVALID(result.outputActions);
      REQUIRE_FALSE(result.outputActions.errors().empty());
      CHECK(result.outputActions.errors()[0].code == testCase.error_code);
      CHECK(result.outputActions.errors()[0].message.find(testCase.diagnostic_fragment) != std::string::npos);
    }
  }
}

TEST_CASE("MTRSim::ComputeODFFilter: Deferred metadata omits unavailable values and shrinks the RGB geometry", "[MTRSim][ComputeODFFilter][Presentation]")
{
  auto data = CreateDeferredMetadataDataStructure();
  auto* namesPtr = data.getDataAs<StringArray>(k_MaterialNamesPath);
  REQUIRE(namesPtr != nullptr);
  namesPtr->setStore(std::make_shared<EmptyStringStore>(ShapeType{2}));
  auto args = ExecutionArguments();
  const auto preflight = ComputeODFFilter{}.preflight(data, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.outputActions);
  for(const auto& value : preflight.outputValues)
  {
    REQUIRE(value.value.find("Deferred until execution") == std::string::npos);
    REQUIRE(value.name != "Section Angles");
    REQUIRE(value.name != "Material");
    REQUIRE(value.name != "Laue Class");
    REQUIRE(value.name != "Euler Plot Limits");
  }
  REQUIRE(GetPreflightValue(preflight, "Selected Phase") == "1");
  REQUIRE(GetPreflightValue(preflight, "Composite Dimensions").find("conservative maximum") != std::string::npos);
  const auto applied = preflight.outputActions.value().applyRegular(data, IDataAction::Mode::Execute);
  SIMPLNX_RESULT_REQUIRE_VALID(applied);
  const auto* geomPtr = data.getDataAs<ImageGeom>(k_OutputGeometryPath);
  REQUIRE(geomPtr != nullptr);
  const auto conservativeDimensions = geomPtr->getDimensions();
  auto* crystalsPtr = data.getDataAs<UInt32Array>(k_CrystalStructuresPath);
  REQUIRE(crystalsPtr != nullptr);
  const auto replaced = crystalsPtr->setDataStore(std::make_shared<DataStore<uint32>>(ShapeType{2}, ShapeType{1}, ebsdlib::CrystalStructure::Hexagonal_High));
  SIMPLNX_RESULT_REQUIRE_VALID(replaced);
  REQUIRE(data.removeData(k_MaterialNamesPath));
  auto* ensemblePtr = data.getDataAs<AttributeMatrix>(k_MaterialNamesPath.getParent());
  REQUIRE(ensemblePtr != nullptr);
  REQUIRE(StringArray::CreateWithValues(data, "Material Names", {2}, {"Invalid", "Titanium (Alpha)"}, ensemblePtr->getId()) != nullptr);
  // This adapter models execution after an upstream reader fills preflight placeholders.
  struct MaterializedComputeFilter : ComputeODFFilter
  {
    using ComputeODFFilter::executeImpl;
  };
  TemporaryOutput output;
  args.insertOrAssign(ComputeODFFilter::k_WritePng_Key, true);
  args.insertOrAssign(ComputeODFFilter::k_OutputPngPath_Key, output.png);
  const auto result = MaterializedComputeFilter{}.executeImpl(data, args, nullptr, {}, false, {});
  SIMPLNX_RESULT_REQUIRE_VALID(result);
  REQUIRE(geomPtr->getNumXCells() == conservativeDimensions[0]);
  REQUIRE(geomPtr->getNumYCells() < conservativeDimensions[1]);
  const auto* rgbPtr = data.getDataAs<UInt8Array>(k_RgbPath);
  REQUIRE(rgbPtr != nullptr);
  REQUIRE(rgbPtr->getTupleShape() == geomPtr->getCellData()->getShape());
  auto decoded = ReadPng(output.png);
  const auto* decodedGeomPtr = decoded.getDataAs<ImageGeom>(DataPath({"Decoded"}));
  const auto* rgbaPtr = decoded.getDataAs<UInt8Array>(DataPath({"Decoded", "Cells", "Pixels"}));
  REQUIRE(decodedGeomPtr != nullptr);
  REQUIRE(rgbaPtr != nullptr);
  REQUIRE(decodedGeomPtr->getDimensions() == geomPtr->getDimensions());
  bool matches = true;
  for(usize tupleIdx = 0; tupleIdx < rgbPtr->getNumberOfTuples(); ++tupleIdx)
  {
    for(usize compIdx = 0; compIdx < 3; ++compIdx)
    {
      matches &= (*rgbPtr)[tupleIdx * 3 + compIdx] == (*rgbaPtr)[tupleIdx * 4 + compIdx];
    }
  }
  REQUIRE(matches);
  UnitTest::CheckArraysInheritTupleDims(data);
}

TEST_CASE("MTRSim::ComputeODFFilter: Placeholder Material Names preserves known crystal previews", "[MTRSim][ComputeODFFilter][Presentation]")
{
  auto data = CreateValidDataStructure();
  const auto args = ExecutionArguments();
  const auto known = ComputeODFFilter{}.preflight(data, args);
  SIMPLNX_RESULT_REQUIRE_VALID(known.outputActions);
  auto* namesPtr = data.getDataAs<StringArray>(k_MaterialNamesPath);
  REQUIRE(namesPtr != nullptr);
  namesPtr->setStore(std::make_shared<EmptyStringStore>(ShapeType{2}));
  IFilter::PreflightResult preflight;
  REQUIRE_NOTHROW(preflight = ComputeODFFilter{}.preflight(data, args));
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.outputActions);
  for(const auto& value : preflight.outputValues)
  {
    REQUIRE(value.name != "Material");
    REQUIRE(value.value.find("Deferred until execution") == std::string::npos);
  }
  REQUIRE(GetPreflightValue(preflight, "Section Angles") == "0, 10, 20, 30, 40, 50 degrees");
  REQUIRE(GetPreflightValue(preflight, "Laue Class") == "Hexagonal 6/mmm (D6h)");
  REQUIRE(GetPreflightValue(preflight, "Euler Plot Limits") == GetPreflightValue(known, "Euler Plot Limits"));
  REQUIRE(GetPreflightValue(preflight, "Composite Dimensions") == GetPreflightValue(known, "Composite Dimensions"));
  REQUIRE(preflight.outputActions.value().actions.size() == known.outputActions.value().actions.size());
}

TEST_CASE("MTRSim::ComputeODFFilter: Placeholder Material Names returns a contextual execution error", "[MTRSim][ComputeODFFilter][Presentation]")
{
  auto data = CreateValidDataStructure();
  auto* namesPtr = data.getDataAs<StringArray>(k_MaterialNamesPath);
  REQUIRE(namesPtr != nullptr);
  namesPtr->setStore(std::make_shared<EmptyStringStore>(ShapeType{2}));
  auto args = ExecutionArguments();
  TemporaryOutput output;
  args.insertOrAssign(ComputeODFFilter::k_WritePng_Key, true);
  args.insertOrAssign(ComputeODFFilter::k_OutputPngPath_Key, output.png);
  IFilter::ExecuteResult result;
  REQUIRE_NOTHROW(result = ComputeODFFilter{}.execute(data, args));
  SIMPLNX_RESULT_REQUIRE_INVALID(result.result);
  REQUIRE(result.result.errors().front().code == -12323);
  const auto& message = result.result.errors().front().message;
  REQUIRE(message.find(k_MaterialNamesPath.toString()) != std::string::npos);
  REQUIRE(message.find("phase 1") != std::string::npos);
  REQUIRE(message.find("not available") != std::string::npos);
  REQUIRE_FALSE(std::filesystem::exists(output.png));
}

TEST_CASE("MTRSim::ComputeODFFilter: Internal presentation adapts Float32 and Float64 identically", "[MTRSim][ComputeODFFilter][Presentation]")
{
  auto data32 = ComputeNumericData();
  auto data64 = ComputeNumericData();
  auto* source32Ptr = data32.getDataAs<Float64Array>(k_OdfArrayPath);
  auto* source64Ptr = data64.getDataAs<Float64Array>(k_OdfArrayPath);
  REQUIRE(source32Ptr != nullptr);
  REQUIRE(source64Ptr != nullptr);
  auto* cellsPtr = data32.getDataAs<AttributeMatrix>(k_OdfArrayPath.getParent());
  REQUIRE(cellsPtr != nullptr);
  const auto floatPath = k_OdfArrayPath.replaceName("Float ODF");
  auto* odf32Ptr = Float32Array::Create(data32, "Float ODF", DataStoreUtilities::CreateDataStore<float32>(data32, floatPath, source32Ptr->getTupleShape(), {1}), cellsPtr->getId());
  REQUIRE(odf32Ptr != nullptr);
  for(usize index = 0; index < odf32Ptr->getSize(); ++index)
  {
    (*odf32Ptr)[index] = static_cast<float32>((*source32Ptr)[index]);
    (*source64Ptr)[index] = (*odf32Ptr)[index];
  }
  auto values = ExecutionInputValues();
  const auto result64 = RenderExisting(data64, values);
  SIMPLNX_RESULT_REQUIRE_VALID(result64);
  values.odfArrayPath = floatPath;
  const auto result32 = RenderExisting(data32, values);
  SIMPLNX_RESULT_REQUIRE_VALID(result32);
  const auto* rgb32Ptr = data32.getDataAs<UInt8Array>(k_RgbPath);
  const auto* rgb64Ptr = data64.getDataAs<UInt8Array>(k_RgbPath);
  REQUIRE(rgb32Ptr != nullptr);
  REQUIRE(rgb64Ptr != nullptr);
  REQUIRE(std::equal(rgb32Ptr->begin(), rgb32Ptr->end(), rgb64Ptr->begin(), rgb64Ptr->end()));
  UnitTest::CheckArraysInheritTupleDims(data32);
  UnitTest::CheckArraysInheritTupleDims(data64);
}

TEST_CASE("MTRSim::ComputeODFFilter: Invalid presentation values preserve an existing PNG destination", "[MTRSim][ComputeODFFilter][Presentation]")
{
  const double value = GENERATE(-2.5, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity());
  TemporaryOutput output;
  const std::string original = "Existing destination must survive a failed render.\n";
  {
    std::ofstream stream(output.png, std::ios::binary);
    stream << original;
    REQUIRE(stream.good());
  }
  auto data = ComputeNumericData();
  auto* odfPtr = data.getDataAs<Float64Array>(k_OdfArrayPath);
  REQUIRE(odfPtr != nullptr);
  // Values outside the rendered crop still require validation.
  (*odfPtr)[odfPtr->getSize() - 1] = value;
  auto values = ExecutionInputValues();
  values.writePng = true;
  values.outputPngPath = output.png;
  const auto result = RenderExisting(data, values);
  SIMPLNX_RESULT_REQUIRE_INVALID(result);
  const auto& message = result.errors().front().message;
  REQUIRE(message.find(k_OdfArrayPath.toString()) != std::string::npos);
  REQUIRE(message.find("phase 1") != std::string::npos);
  REQUIRE(message.find("index 186623") != std::string::npos);
  std::ifstream stream(output.png, std::ios::binary);
  REQUIRE(std::string(std::istreambuf_iterator<char>(stream), {}) == original);
  REQUIRE(std::distance(std::filesystem::directory_iterator(output.directory), std::filesystem::directory_iterator{}) == 1);
}

TEST_CASE("MTRSim::ComputeODFFilter: Empty selected phase rejects presentation", "[MTRSim][ComputeODFFilter][Presentation]")
{
  UnitTest::LoadPlugins();
  const bool masked = GENERATE(false, true);
  const auto scale = GENERATE(ChoicesParameter::ValueType{0}, ChoicesParameter::ValueType{1});
  const auto mode = GENERATE(std::pair{false, false}, std::pair{false, true}, std::pair{true, false}, std::pair{true, true});
  const auto units = GENERATE(ChoicesParameter::ValueType{0}, ChoicesParameter::ValueType{1});
  CAPTURE(masked, scale, mode.first, mode.second, units);
  auto data = CreateValidDataStructure();
  auto* phasesPtr = data.getDataAs<Int32Array>(DataPath({"Input Cells", "Phases"}));
  REQUIRE(phasesPtr != nullptr);
  auto args = ExecutionArguments();
  if(masked)
  {
    const DataPath maskPath({"Input Cells", "Mask"});
    auto* maskPtr = BoolArray::Create(data, "Mask", DataStoreUtilities::CreateDataStore<bool>(data, maskPath, {8}, {1}), phasesPtr->getParentIds().front());
    REQUIRE(maskPtr != nullptr);
    std::fill(maskPtr->begin(), maskPtr->end(), false);
    args.insertOrAssign(ComputeODFFilter::k_UseMask_Key, true);
    args.insertOrAssign(ComputeODFFilter::k_MaskArrayPath_Key, maskPath);
  }
  else
  {
    std::fill(phasesPtr->begin(), phasesPtr->end(), 0);
  }
  TemporaryOutput output;
  const std::string original = "Existing PNG destination";
  {
    std::ofstream stream(output.png, std::ios::binary);
    stream << original;
    REQUIRE(stream.good());
  }
  args.insertOrAssign(ComputeODFFilter::k_ScaleMode_Key, scale);
  args.insertOrAssign(ComputeODFFilter::k_OutputUnits_Key, units);
  args.insertOrAssign(ComputeODFFilter::k_ManualMinimum_Key, 0.25);
  args.insertOrAssign(ComputeODFFilter::k_ManualMaximum_Key, 2.5);
  args.insertOrAssign(ComputeODFFilter::k_WritePng_Key, mode.first);
  args.insertOrAssign(ComputeODFFilter::k_SaveAsImageGeometry_Key, mode.second);
  args.insertOrAssign(ComputeODFFilter::k_OutputPngPath_Key, output.png);
  const auto result = ComputeODFFilter{}.execute(data, args);
  if(mode.first || mode.second)
  {
    CHECK(result.result.invalid());
    if(result.result.invalid())
    {
      CHECK(result.result.errors().front().message.find("phase 1") != std::string::npos);
      CHECK(result.result.errors().front().message.find("no voxels contributed") != std::string::npos);
    }
    const auto* rgbPtr = data.getDataAs<UInt8Array>(k_RgbPath);
    CHECK((rgbPtr == nullptr || std::all_of(rgbPtr->begin(), rgbPtr->end(), [](uint8 value) { return value == 0; })));
  }
  else
  {
    SIMPLNX_RESULT_REQUIRE_VALID(result.result);
    const auto* odfPtr = data.getDataAs<Float64Array>(k_OdfArrayPath);
    REQUIRE(odfPtr != nullptr);
    CHECK(std::all_of(odfPtr->begin(), odfPtr->end(), [](float64 value) { return value == 0.0; }));
    TemporaryOutput persistence;
    const auto numericFile = persistence.directory / "empty_phase.dream3d";
    UnitTest::WriteTestDataStructure(data, numericFile);
    const auto restored = UnitTest::LoadDataStructure(numericFile);
    const auto* metadataPtr = restored.getDataAs<StringArray>(DataPath({"ODF", "ODF Metadata", "ODF Values Units"}));
    REQUIRE(metadataPtr != nullptr);
    CHECK((*metadataPtr)[0] == (units == 0 ? "Count-Density" : "MUD"));
    CHECK(data.getDataAs<UInt8Array>(k_RgbPath) == nullptr);
  }
  std::ifstream stream(output.png, std::ios::binary);
  CHECK(std::string(std::istreambuf_iterator<char>(stream), {}) == original);
  CHECK(std::distance(std::filesystem::directory_iterator(output.directory), std::filesystem::directory_iterator{}) == 1);
  UnitTest::CheckArraysInheritTupleDims(data);
}

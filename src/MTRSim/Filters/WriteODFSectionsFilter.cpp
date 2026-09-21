#include "WriteODFSectionsFilter.hpp"

#include "MTRSim/Filters/Algorithms/WriteODFSections.hpp"

#include "EbsdLib/LaueOps/LaueOps.h"
#include "EbsdLib/Utilities/ODFSectionCompositor.h"

#include "simplnx/Common/TypesUtility.hpp"
#include "simplnx/DataStructure/Geometry/IGeometry.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/IDataArray.hpp"
#include "simplnx/DataStructure/IDataStore.hpp"
#include "simplnx/DataStructure/StringArray.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/CreateColorMapParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataPathSelectionParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"
#include "simplnx/Utilities/ColorTableUtilities.hpp"
#include "simplnx/Utilities/StringUtilities.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <cmath>
#include <vector>

using namespace nx::core;

namespace
{
constexpr StringLiteral k_CellDataName = "Cell Data";
constexpr StringLiteral k_RgbArrayName = "ODF Sections RGB";
} // namespace

std::string WriteODFSectionsFilter::name() const
{
  return FilterTraits<WriteODFSectionsFilter>::name.str();
}

std::string WriteODFSectionsFilter::className() const
{
  return FilterTraits<WriteODFSectionsFilter>::className;
}

Uuid WriteODFSectionsFilter::uuid() const
{
  return FilterTraits<WriteODFSectionsFilter>::uuid;
}

std::string WriteODFSectionsFilter::humanName() const
{
  return "Write ODF Section Image";
}

std::vector<std::string> WriteODFSectionsFilter::defaultTags() const
{
  return {className(), "MTRSim", "ODF", "Sections", "Image", "PNG"};
}

Parameters WriteODFSectionsFilter::parameters() const
{
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input ODF"});
  params.insert(std::make_unique<GeometrySelectionParameter>(k_InputImageGeometryPath_Key, "ODF Image Geometry", "Image geometry that contains the full Bunge Euler grid.", DataPath{},
                                                             GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_ODFArrayPath_Key, "ODF Array", "Scalar Float32 or Float64 ODF values on the selected geometry cell data.", DataPath{},
                                                          ArraySelectionParameter::AllowedTypes{DataType::float32, DataType::float64}, ArraySelectionParameter::AllowedComponentShapes{{1}}));
  params.insert(std::make_unique<ChoicesParameter>(k_InputUnits_Key, "Input Units", "Units of the selected ODF values.", 0ULL, ChoicesParameter::Choices{"MUD", "Count-Density"}));

  params.insertSeparator(Parameters::Separator{"Phase Metadata"});
  params.insert(std::make_unique<ArraySelectionParameter>(k_CrystalStructuresPath_Key, "Crystal Structures", "UInt32 crystal-structure code for each ensemble tuple.", DataPath{},
                                                          ArraySelectionParameter::AllowedTypes{DataType::uint32}, ArraySelectionParameter::AllowedComponentShapes{{1}}));
  params.insert(std::make_unique<DataPathSelectionParameter>(k_MaterialNamesPath_Key, "Material Names", "StringArray material name for each ensemble tuple.", DataPath{}));
  params.insert(std::make_unique<Int32Parameter>(k_PhaseIndex_Key, "Phase Index", "One-based ensemble tuple to render.", 1));

  params.insertSeparator(Parameters::Separator{"Presentation"});
  params.insert(std::make_unique<StringParameter>(k_Title_Key, "Figure Title", "Title drawn above the ODF section montage.", "ODF Sections"));
  params.insert(std::make_unique<Int32Parameter>(k_NumSections_Key, "Number of Sections", "Number of constant-phi2 sections.", 6));
  params.insert(std::make_unique<Int32Parameter>(k_SectionsPerRow_Key, "Sections Per Row", "Number of section panels in each row.", 3));
  params.insert(std::make_unique<Int32Parameter>(k_SectionWidth_Key, "Section Width", "Width of each ODF section panel in pixels.", 512));
  params.insert(std::make_unique<CreateColorMapParameter>(k_ColorTablePreset_Key, "Color Table Preset", "RGB preset for the shared MUD color scale.", "Rainbow"));
  params.insertLinkableParameter(
      std::make_unique<ChoicesParameter>(k_ScaleMode_Key, "Scale Mode", "Selects an automatic or manual shared MUD range.", 0ULL, ChoicesParameter::Choices{"Automatic", "Manual"}));
  params.insert(std::make_unique<Float64Parameter>(k_ManualMinimum_Key, "Manual Minimum MUD", "Minimum value for manual MUD scaling.", 0.0));
  params.insert(std::make_unique<Float64Parameter>(k_ManualMaximum_Key, "Manual Maximum MUD", "Maximum value for manual MUD scaling.", 1.0));
  params.linkParameters(k_ScaleMode_Key, k_ManualMinimum_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  params.linkParameters(k_ScaleMode_Key, k_ManualMaximum_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));

  params.insertSeparator(Parameters::Separator{"Output"});
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_WritePng_Key, "Write PNG", "Writes the annotated montage to a PNG file.", true));
  params.insert(std::make_unique<FileSystemPathParameter>(k_OutputPngPath_Key, "Output PNG", "Output file for the annotated montage.", std::filesystem::path{},
                                                          FileSystemPathParameter::ExtensionsType{".png"}, FileSystemPathParameter::PathType::OutputFile));
  params.linkParameters(k_WritePng_Key, k_OutputPngPath_Key, true);
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_SaveAsImageGeometry_Key, "Save as Image Geometry", "Creates an RGB Image Geometry for the montage.", true));
  params.insert(std::make_unique<DataGroupCreationParameter>(k_OutputImageGeometryPath_Key, "Output Image Geometry", "Path for the RGB montage geometry.", DataPath({"ODF Sections"})));
  params.linkParameters(k_SaveAsImageGeometry_Key, k_OutputImageGeometryPath_Key, true);
  return params;
}

IFilter::VersionType WriteODFSectionsFilter::parametersVersion() const
{
  return 1;
}

IFilter::UniquePointer WriteODFSectionsFilter::clone() const
{
  return std::make_unique<WriteODFSectionsFilter>();
}

IFilter::PreflightResult WriteODFSectionsFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                               const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  const auto imageGeometryPath = filterArgs.value<DataPath>(k_InputImageGeometryPath_Key);
  const auto odfArrayPath = filterArgs.value<DataPath>(k_ODFArrayPath_Key);
  const auto inputUnits = filterArgs.value<ChoicesParameter::ValueType>(k_InputUnits_Key);
  const auto crystalStructuresPath = filterArgs.value<DataPath>(k_CrystalStructuresPath_Key);
  const auto materialNamesPath = filterArgs.value<DataPath>(k_MaterialNamesPath_Key);
  const int32 phaseIndex = filterArgs.value<int32>(k_PhaseIndex_Key);
  const int32 sectionCount = filterArgs.value<int32>(k_NumSections_Key);
  const int32 sectionsPerRow = filterArgs.value<int32>(k_SectionsPerRow_Key);
  const int32 sectionWidth = filterArgs.value<int32>(k_SectionWidth_Key);
  const auto colorTablePreset = filterArgs.value<CreateColorMapParameter::ValueType>(k_ColorTablePreset_Key);
  const auto scaleMode = filterArgs.value<ChoicesParameter::ValueType>(k_ScaleMode_Key);
  const float64 manualMinimum = filterArgs.value<float64>(k_ManualMinimum_Key);
  const float64 manualMaximum = filterArgs.value<float64>(k_ManualMaximum_Key);
  const bool writePng = filterArgs.value<bool>(k_WritePng_Key);
  const auto outputPngPath = filterArgs.value<FileSystemPathParameter::ValueType>(k_OutputPngPath_Key);
  const bool saveAsImageGeometry = filterArgs.value<bool>(k_SaveAsImageGeometry_Key);
  const auto outputImageGeometryPath = filterArgs.value<DataPath>(k_OutputImageGeometryPath_Key);

  if(!writePng && !saveAsImageGeometry)
  {
    return {MakeErrorResult<OutputActions>(-12300, "At least one output must be enabled. Write PNG=false and Save as Image Geometry=false.")};
  }

  const auto& imageGeom = dataStructure.getDataRefAs<ImageGeom>(imageGeometryPath);
  const auto& odfArray = dataStructure.getDataRefAs<IDataArray>(odfArrayPath);
  const auto& crystalStructures = dataStructure.getDataRefAs<UInt32Array>(crystalStructuresPath);
  const auto* materialNames = dataStructure.getDataAs<StringArray>(materialNamesPath);

  if((odfArray.getDataType() != DataType::float32 && odfArray.getDataType() != DataType::float64) || odfArray.getNumberOfComponents() != 1)
  {
    return {MakeErrorResult<OutputActions>(-12301, fmt::format("ODF array '{}' must be a scalar Float32 or Float64 array; selected type '{}' with {} components.", odfArrayPath.toString(),
                                                               DataTypeToString(odfArray.getDataType()), odfArray.getNumberOfComponents()))};
  }

  if(imageGeom.getCellData() == nullptr || odfArrayPath.getParent() != imageGeom.getCellDataPath())
  {
    const std::string cellDataPath = imageGeom.getCellData() == nullptr ? "(none)" : imageGeom.getCellDataPath().toString();
    return {MakeErrorResult<OutputActions>(-12302, fmt::format("ODF array '{}' must be directly under the selected Image Geometry '{}' cell Attribute Matrix '{}'.", odfArrayPath.toString(),
                                                               imageGeometryPath.toString(), cellDataPath))};
  }

  const auto origin = imageGeom.getOrigin();
  if(origin[0] != 0.0F || origin[1] != 0.0F || origin[2] != 0.0F)
  {
    return {
        MakeErrorResult<OutputActions>(-12303, fmt::format("ODF Image Geometry '{}' origin ({}, {}, {}) degrees must be (0, 0, 0).", imageGeometryPath.toString(), origin[0], origin[1], origin[2]))};
  }

  const auto spacing = imageGeom.getSpacing();
  if(!std::isfinite(spacing[0]) || !std::isfinite(spacing[1]) || !std::isfinite(spacing[2]) || spacing[0] <= 0.0F || spacing[1] <= 0.0F || spacing[2] <= 0.0F || spacing[0] != spacing[1] ||
     spacing[0] != spacing[2])
  {
    return {MakeErrorResult<OutputActions>(
        -12304, fmt::format("ODF Image Geometry '{}' spacing ({}, {}, {}) degrees must be finite, positive, and uniform.", imageGeometryPath.toString(), spacing[0], spacing[1], spacing[2]))};
  }

  constexpr float64 k_BinCountTolerance = 1.0e-6;
  const bool coversFullBungeCube = imageGeom.getNumXCells() > 0 && imageGeom.getNumYCells() > 0 && imageGeom.getNumZCells() > 0 &&
                                   std::abs(360.0 / spacing[0] - static_cast<float64>(imageGeom.getNumXCells())) <= k_BinCountTolerance &&
                                   std::abs(180.0 / spacing[1] - static_cast<float64>(imageGeom.getNumYCells())) <= k_BinCountTolerance &&
                                   std::abs(360.0 / spacing[2] - static_cast<float64>(imageGeom.getNumZCells())) <= k_BinCountTolerance;
  if(!coversFullBungeCube)
  {
    return {MakeErrorResult<OutputActions>(-12305, fmt::format("ODF Image Geometry '{}' dimensions ({}) with spacing ({}, {}, {}) degrees must cover the full 360 by 180 by 360 degree Bunge cube.",
                                                               imageGeometryPath.toString(), StringUtilities::formatDimensions3D(imageGeom.getDimensions()), spacing[0], spacing[1], spacing[2]))};
  }

  if(odfArray.getNumberOfTuples() != imageGeom.getNumberOfCells())
  {
    return {MakeErrorResult<OutputActions>(-12301, fmt::format("ODF array '{}' contains {} tuples, but Image Geometry '{}' contains {} cells.", odfArrayPath.toString(), odfArray.getNumberOfTuples(),
                                                               imageGeometryPath.toString(), imageGeom.getNumberOfCells()))};
  }

  if(materialNames == nullptr)
  {
    return {MakeErrorResult<OutputActions>(-12306, fmt::format("Material Names path '{}' must select a StringArray; the selected object has a different type.", materialNamesPath.toString()))};
  }
  if(crystalStructures.getNumberOfTuples() != materialNames->getNumberOfTuples())
  {
    return {
        MakeErrorResult<OutputActions>(-12306, fmt::format("Crystal Structures '{}' contains {} tuples, but Material Names '{}' contains {} tuples. The ensemble tuple counts must match.",
                                                           crystalStructuresPath.toString(), crystalStructures.getNumberOfTuples(), materialNamesPath.toString(), materialNames->getNumberOfTuples()))};
  }
  const usize ensembleTupleCount = crystalStructures.getNumberOfTuples();
  if(phaseIndex <= 0 || static_cast<usize>(phaseIndex) >= ensembleTupleCount)
  {
    return {MakeErrorResult<OutputActions>(-12307, fmt::format("Phase Index ({}) must be in [1, {}] for Crystal Structures '{}' and Material Names '{}'.", phaseIndex,
                                                               ensembleTupleCount == 0 ? 0 : ensembleTupleCount - 1, crystalStructuresPath.toString(), materialNamesPath.toString()))};
  }

  const auto metadataStoreType = crystalStructures.getIDataStoreRef().getStoreType();
  const bool deferredMetadata = metadataStoreType == IDataStore::StoreType::Empty;
  // Upstream preflight actions provide tuple shapes but no phase values. Triclinic bounds give the largest page.
  const uint32 crystalStructure = deferredMetadata ? ebsdlib::CrystalStructure::Triclinic : crystalStructures[static_cast<usize>(phaseIndex)];
  const auto limits = ebsdlib::GetODFEulerPlotLimits(crystalStructure);
  if(limits.phi1MaxDeg <= 0.0 || limits.phiMaxDeg <= 0.0 || limits.phi2MaxDeg <= 0.0)
  {
    return {MakeErrorResult<OutputActions>(
        -12308, fmt::format("Crystal Structures '{}' contains unsupported crystal-structure value {} at phase index {}.", crystalStructuresPath.toString(), crystalStructure, phaseIndex))};
  }
  if(sectionCount < 2)
  {
    return {MakeErrorResult<OutputActions>(-12309, fmt::format("Number of Sections ({}) must be at least 2.", sectionCount))};
  }
  if(sectionsPerRow < 1 || sectionsPerRow > sectionCount)
  {
    return {MakeErrorResult<OutputActions>(-12310, fmt::format("Sections Per Row ({}) must be in [1, Number of Sections ({})].", sectionsPerRow, sectionCount))};
  }
  if(sectionWidth < 128)
  {
    return {MakeErrorResult<OutputActions>(-12311, fmt::format("Section Width ({}) pixels must be at least 128 pixels.", sectionWidth))};
  }

  auto colorControls = ColorTableUtilities::ExtractControlPoints(colorTablePreset);
  if(colorControls.invalid())
  {
    return {MakeErrorResult<OutputActions>(-12312, fmt::format("Color Table Preset '{}' is not a valid RGB preset.", colorTablePreset))};
  }
  if(scaleMode == 1 && (!std::isfinite(manualMinimum) || !std::isfinite(manualMaximum) || manualMinimum < 0.0 || manualMinimum >= manualMaximum))
  {
    return {MakeErrorResult<OutputActions>(-12313, fmt::format("Manual MUD range must have finite endpoints and satisfy 0 <= minimum < maximum; received [{}, {}].", manualMinimum, manualMaximum))};
  }
  if(writePng && outputPngPath.extension() != ".png")
  {
    return {MakeErrorResult<OutputActions>(-12314, fmt::format("Output PNG file '{}' has extension '{}'; use the '.png' extension.", outputPngPath.string(), outputPngPath.extension().string()))};
  }

  ebsdlib::ODFSectionConfiguration config;
  // Image geometry XYZ corresponds to phi2, PHI, phi1; the borrowed values retain phi2-fastest storage.
  config.grid.dimensions = {imageGeom.getNumZCells(), imageGeom.getNumYCells(), imageGeom.getNumXCells()};
  config.grid.originDeg = {origin[2], origin[1], origin[0]};
  config.grid.spacingDeg = {spacing[2], spacing[1], spacing[0]};
  config.grid.laueOpsIndex = crystalStructure;
  config.sectionCount = static_cast<usize>(sectionCount);
  config.sectionsPerRow = static_cast<usize>(sectionsPerRow);
  config.sectionWidth = sectionWidth;
  config.scaleMode = scaleMode == 0 ? ebsdlib::ODFScaleMode::Automatic : ebsdlib::ODFScaleMode::Manual;
  config.manualMinimumMUD = manualMinimum;
  config.manualMaximumMUD = manualMaximum;
  config.colorControlPoints = std::move(colorControls.value());
  config.title = filterArgs.value<StringParameter::ValueType>(k_Title_Key);
  if(!deferredMetadata)
  {
    config.materialName = (*materialNames)[static_cast<usize>(phaseIndex)];
  }
  config.phaseNumber = phaseIndex;

  const auto layout = ebsdlib::ComputeODFSectionLayout(config);
  if(layout.pageWidth <= 0 || layout.pageHeight <= 0)
  {
    return {MakeErrorResult<OutputActions>(
        -12315, fmt::format("ODF section layout is not representable for Section Width {}, Number of Sections {}, and Sections Per Row {}.", sectionWidth, sectionCount, sectionsPerRow))};
  }

  Result<OutputActions> actions;
  if(saveAsImageGeometry)
  {
    actions.value().appendAction(
        std::make_unique<CreateImageGeometryAction>(outputImageGeometryPath, CreateImageGeometryAction::DimensionType{static_cast<usize>(layout.pageWidth), static_cast<usize>(layout.pageHeight), 1},
                                                    CreateImageGeometryAction::OriginType{0.0F, 0.0F, 0.0F}, CreateImageGeometryAction::SpacingType{1.0F, 1.0F, 1.0F}, k_CellDataName.str()));
    const DataPath rgbArrayPath = outputImageGeometryPath.createChildPath(k_CellDataName.str()).createChildPath(k_RgbArrayName.str());
    actions.value().appendAction(
        std::make_unique<CreateArrayAction>(DataType::uint8, ShapeType{1, static_cast<usize>(layout.pageHeight), static_cast<usize>(layout.pageWidth)}, ShapeType{3}, rgbArrayPath));
  }

  const auto laueNames = ebsdlib::LaueOps::GetLaueNames();
  const auto sectionAngles = deferredMetadata ? std::vector<double>{} : ebsdlib::GenerateODFSectionAngles(crystalStructure, config.sectionCount);
  std::vector<PreflightValue> values;
  values.push_back(
      {"Section Angles", deferredMetadata ? std::string{"Deferred until execution; phase metadata is not available during preflight."} : fmt::format("{} degrees", fmt::join(sectionAngles, ", "))});
  values.push_back({"Phase", std::to_string(phaseIndex)});
  if(deferredMetadata)
  {
    values.push_back({"Selected Phase Metadata", fmt::format("Deferred until execution for phase {} from '{}' and '{}'.", phaseIndex, crystalStructuresPath.toString(), materialNamesPath.toString())});
    values.push_back({"Laue Group", "Deferred until execution; the Triclinic envelope is used only for maximum layout sizing."});
  }
  else
  {
    values.push_back({"Material", config.materialName});
    values.push_back({"Laue Class", laueNames[crystalStructure]});
  }
  values.push_back({"Input Units", inputUnits == 0 ? "MUD" : "Count-Density"});
  values.push_back({"Rendered Units", "MUD"});
  values.push_back({"Full Grid Dimensions (phi1 x PHI x phi2)", fmt::format("{} x {} x {}", config.grid.dimensions[0], config.grid.dimensions[1], config.grid.dimensions[2])});
  values.push_back({"Full Grid Spacing (phi1, PHI, phi2)", fmt::format("{}, {}, {} degrees", config.grid.spacingDeg[0], config.grid.spacingDeg[1], config.grid.spacingDeg[2])});
  values.push_back({"Euler Plot Limits", deferredMetadata ? std::string{"Deferred until execution; maximum envelope is 360 x 180 x 360 degrees."} :
                                                            fmt::format("phi1 < {}, PHI < {}, phi2 < {} degrees", limits.phi1MaxDeg, limits.phiMaxDeg, limits.phi2MaxDeg)});
  values.push_back({"Section Layout", fmt::format("{} columns x {} rows", layout.columns, layout.rows)});
  values.push_back({"Composite Dimensions", deferredMetadata ?
                                                fmt::format("{} x {} pixels (conservative maximum estimate; exact dimensions are set during execution).", layout.pageWidth, layout.pageHeight) :
                                                fmt::format("{} x {} pixels", layout.pageWidth, layout.pageHeight)});
  values.push_back({"Scale Mode", scaleMode == 0 ? std::string{"Automatic"} : fmt::format("Manual: [{}, {}] MUD", manualMinimum, manualMaximum)});
  return {std::move(actions), std::move(values)};
}

Result<> WriteODFSectionsFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                             const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  WriteODFSectionsInputValues inputValues;
  inputValues.inputImageGeometryPath = filterArgs.value<DataPath>(k_InputImageGeometryPath_Key);
  inputValues.odfArrayPath = filterArgs.value<DataPath>(k_ODFArrayPath_Key);
  inputValues.inputUnits = filterArgs.value<ChoicesParameter::ValueType>(k_InputUnits_Key);
  inputValues.crystalStructuresPath = filterArgs.value<DataPath>(k_CrystalStructuresPath_Key);
  inputValues.materialNamesPath = filterArgs.value<DataPath>(k_MaterialNamesPath_Key);
  inputValues.phaseIndex = filterArgs.value<int32>(k_PhaseIndex_Key);
  inputValues.title = filterArgs.value<StringParameter::ValueType>(k_Title_Key);
  inputValues.numberOfSections = filterArgs.value<int32>(k_NumSections_Key);
  inputValues.sectionsPerRow = filterArgs.value<int32>(k_SectionsPerRow_Key);
  inputValues.sectionWidth = filterArgs.value<int32>(k_SectionWidth_Key);
  inputValues.colorTablePreset = filterArgs.value<CreateColorMapParameter::ValueType>(k_ColorTablePreset_Key);
  inputValues.scaleMode = filterArgs.value<ChoicesParameter::ValueType>(k_ScaleMode_Key);
  inputValues.manualMinimumMUD = filterArgs.value<float64>(k_ManualMinimum_Key);
  inputValues.manualMaximumMUD = filterArgs.value<float64>(k_ManualMaximum_Key);
  inputValues.writePng = filterArgs.value<bool>(k_WritePng_Key);
  inputValues.outputPngPath = filterArgs.value<FileSystemPathParameter::ValueType>(k_OutputPngPath_Key);
  inputValues.saveAsImageGeometry = filterArgs.value<bool>(k_SaveAsImageGeometry_Key);
  inputValues.outputImageGeometryPath = filterArgs.value<DataPath>(k_OutputImageGeometryPath_Key);
  return WriteODFSections(dataStructure, messageHandler, shouldCancel, &inputValues)();
}

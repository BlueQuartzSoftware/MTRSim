#include "ComputeODFFilter.hpp"

#include "MTRSim/Filters/Algorithms/ComputeODF.hpp"
#include "MTRSim/Filters/Algorithms/GenerateODFSections.hpp"
#include "MTRSim/Filters/ODFMetadata.hpp"

#include "EbsdLib/Utilities/ODFSectionCompositor.h"
#include "fmt/ranges.h"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataGroup.hpp"
#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/Geometry/IGeometry.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/IDataArray.hpp"
#include "simplnx/DataStructure/StringArray.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateDataGroupAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Filter/Actions/CreateStringArrayAction.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/CreateColorMapParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/DataPathSelectionParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"
#include "simplnx/Utilities/ColorTableUtilities.hpp"

#include <EbsdLib/LaueOps/LaueOps.h>

#include <fmt/format.h>

#include <cmath>
#include <cstdint>
#include <vector>

using namespace nx::core;

namespace nx::core
{
//------------------------------------------------------------------------------
std::string ComputeODFFilter::name() const
{
  return FilterTraits<ComputeODFFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string ComputeODFFilter::className() const
{
  return FilterTraits<ComputeODFFilter>::className;
}

//------------------------------------------------------------------------------
Uuid ComputeODFFilter::uuid() const
{
  return FilterTraits<ComputeODFFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string ComputeODFFilter::humanName() const
{
  return "Compute ODF From Euler Angles";
}

//------------------------------------------------------------------------------
std::vector<std::string> ComputeODFFilter::defaultTags() const
{
  return {className(), "MTRSim", "ODF", "EBSD", "Compute", "Crystallography", "Statistics"};
}

//------------------------------------------------------------------------------
Parameters ComputeODFFilter::parameters() const
{
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<Int32Parameter>(k_PhaseIndex_Key, "Phase Index", "Positive ensemble index of the phase to include in the ODF. All other phase labels are excluded.", 1));
  params.insert(std::make_unique<BoolParameter>(k_ApplySmoothing_Key, "Apply Smoothing",
                                                "When true, each symmetric Euler tuple is distributed across 27 bins (1 "
                                                "center + 6 faces + 12 edges + 8 corners) using "
                                                "the MATLAB calc_ODF.m tri-linear weights. When false, each tuple "
                                                "deposits 1.0 into its center bin only.",
                                                true));
  params.insert(std::make_unique<Float32Parameter>(k_BinSizeDeg_Key, "Bin Size [deg]",
                                                   "Uniform bin size in degrees, used for all three Bunge-Euler axes (phi1, "
                                                   "PHI, phi2). 360/binSize and 180/binSize "
                                                   "must both yield a positive integer (e.g. 5.0 -> 72 x 36 x 72 bins).",
                                                   5.0f));
  params.insert(std::make_unique<ChoicesParameter>(k_OutputUnits_Key, "Output Units",
                                                   "Count-Density emits the raw normalized histogram (sums to 1 over all "
                                                   "Bunge bins; bit-exact match to MATLAB calc_ODF.m). "
                                                   "MUD applies the per-PHI-row sin(PHI) Jacobian correction so values "
                                                   "represent density on SO(3) divided by uniform "
                                                   "density (matches MTEX plot(odf) convention; sub-uniform regions read < "
                                                   "1, peak texture reads several MUD).",
                                                   0ULL, ChoicesParameter::Choices{"Count-Density", "MUD"}));

  params.insertSeparator(Parameters::Separator{"Required Input Cell Data"});
  params.insert(std::make_unique<ArraySelectionParameter>(k_EulerAngles_Key, "Cell Euler Angles",
                                                          "Per-voxel Bunge Euler angles (phi1, PHI, phi2) in radians. Must be a "
                                                          "3-component Float32 array on a cell-level "
                                                          "AttributeMatrix.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::float32}, ArraySelectionParameter::AllowedComponentShapes{{3}}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_Phases_Key, "Cell Phases",
                                                          "Per-voxel phase label (1-based; 0 marks excluded voxels). Must be a "
                                                          "1-component Int32 array on the same "
                                                          "AttributeMatrix as the Euler angles.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::int32}, ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insertSeparator(Parameters::Separator{"Optional Data Mask"});
  params.insertLinkableParameter(
      std::make_unique<BoolParameter>(k_UseMask_Key, "Use Mask Array", "When enabled, only voxels of the selected phase where the mask is true contribute to the ODF.", false));
  params.insert(std::make_unique<ArraySelectionParameter>(k_MaskArrayPath_Key, "Mask Array",
                                                          "Per-voxel boolean mask. Must be a 1-component Bool array on the same "
                                                          "AttributeMatrix as the Euler angles. Only "
                                                          "consulted when 'Use Mask Array' is enabled.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::boolean, DataType::uint8}, ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insertSeparator(Parameters::Separator{"Required Input Ensemble Data"});
  params.insert(std::make_unique<DataPathSelectionParameter>(k_MaterialNames_Key, "Material Names", "StringArray with one material name per ensemble tuple, indexed by phase label.", DataPath{}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_CrystalStructures_Key, "Crystal Structures",
                                                          "Per-phase EbsdLib crystal-structure code (e.g. 0 = Hexagonal_High, 1 = "
                                                          "Cubic_High). Must be a 1-component "
                                                          "UInt32 array on the ensemble-level AttributeMatrix; indexed by phase "
                                                          "label.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::uint32}, ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insertSeparator(Parameters::Separator{"Output Mode"});
  params.insertLinkableParameter(std::make_unique<ChoicesParameter>(k_OutputMode_Key, "Output Mode",
                                                                    "Whether to create a new ODF ImageGeometry or append a new component to "
                                                                    "an existing one. In Append mode the bin "
                                                                    "size is derived from the existing geometry's spacing and the user's bin "
                                                                    "size input is ignored.",
                                                                    0ULL, ChoicesParameter::Choices{"Create New ODF Geometry", "Append to Existing ODF Geometry"}));

  params.insertSeparator(Parameters::Separator{"Output Data Object(s)"});
  params.insert(std::make_unique<DataGroupCreationParameter>(k_OutputImageGeometry_Key, "Output ODF Image Geometry",
                                                             "Path at which the ODF ImageGeom will be created. Its (X, Y, Z) "
                                                             "dimensions map to (phi2, PHI, phi1).",
                                                             DataPath({"ODF"})));
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellAttrMatName_Key, "Cell Attribute Matrix Name",
                                                          "Name of the cell AttributeMatrix created under the output ImageGeom. "
                                                          "The Float64 ODF array is placed inside it.",
                                                          "Cell Data"));
  params.insert(std::make_unique<GeometrySelectionParameter>(k_ExistingOdfGeometry_Key, "Existing ODF Image Geometry",
                                                             "ImageGeom representing an existing Euler-space ODF grid that this "
                                                             "filter will append a new component to. Spacing must be "
                                                             "uniform across all three axes.",
                                                             DataPath{}, GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<DataObjectNameParameter>(k_ComponentName_Key, "ODF Component Array Name",
                                                          "Name of the Float64 single-component ODF array created on the output "
                                                          "cell AttributeMatrix.",
                                                          "Component 1"));

  params.linkParameters(k_UseMask_Key, k_MaskArrayPath_Key, true);
  params.linkParameters(k_OutputMode_Key, k_OutputImageGeometry_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  params.linkParameters(k_OutputMode_Key, k_CellAttrMatName_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  params.linkParameters(k_OutputMode_Key, k_BinSizeDeg_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  params.linkParameters(k_OutputMode_Key, k_ExistingOdfGeometry_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));

  params.insertSeparator(Parameters::Separator{"Output Image Generation"});
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_WritePng_Key, "Write PNG", "Writes the ODF to a PNG file.", false));
  params.insert(std::make_unique<FileSystemPathParameter>(k_OutputPngPath_Key, "Output PNG", "Output file for the ODF image.", std::filesystem::path{}, FileSystemPathParameter::ExtensionsType{".png"},
                                                          FileSystemPathParameter::PathType::OutputFile));
  params.linkParameters(k_WritePng_Key, k_OutputPngPath_Key, true);

  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_SaveAsImageGeometry_Key, "Save as Image Geometry", "Creates an RGB Image Geometry for the montage.", true));
  params.insert(std::make_unique<DataGroupCreationParameter>(k_OutputPresentationGeometry_Key, "Output Image Geometry", "Path for the RGB montage geometry.", DataPath({"ODF Sections"})));

  params.insertSeparator(Parameters::Separator{"Presentation"});
  params.insert(std::make_unique<StringParameter>(k_ODF_Image_Title_Key, "ODF Figure Title", "Title drawn above the ODF section montage.", "ODF Sections"));
  params.insert(std::make_unique<Int32Parameter>(k_NumSections_Key, "Number of Sections", "Number of constant-phi2 sections.", 6));
  params.insert(std::make_unique<Int32Parameter>(k_SectionsPerRow_Key, "Sections Per Row", "Number of section panels in each row.", 2));
  params.insert(std::make_unique<Int32Parameter>(k_SectionWidth_Key, "Section Width (pixels)", "Width of each ODF section panel in pixels.", 768));
  params.insert(std::make_unique<CreateColorMapParameter>(k_ColorTablePreset_Key, "Color Table Preset", "RGB preset for the shared MUD color scale.", "Rainbow Blended White"));
  params.insertLinkableParameter(
      std::make_unique<ChoicesParameter>(k_ScaleMode_Key, "Scale Mode", "Selects an automatic or manual shared MUD range.", 0ULL, ChoicesParameter::Choices{"Automatic", "Manual"}));
  params.insert(std::make_unique<Float64Parameter>(k_ManualMinimum_Key, "Manual Minimum MUD", "Minimum value for manual MUD scaling.", 0.0));
  params.insert(std::make_unique<Float64Parameter>(k_ManualMaximum_Key, "Manual Maximum MUD", "Maximum value for manual MUD scaling.", 1.0));
  params.linkParameters(k_ScaleMode_Key, k_ManualMinimum_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));
  params.linkParameters(k_ScaleMode_Key, k_ManualMaximum_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));

  params.linkParameters(k_SaveAsImageGeometry_Key, k_OutputPresentationGeometry_Key, true);

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ComputeODFFilter::parametersVersion() const
{
  return 6;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ComputeODFFilter::clone() const
{
  return std::make_unique<ComputeODFFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ComputeODFFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& /*messageHandler*/,
                                                         const std::atomic_bool& /*shouldCancel*/, const ExecutionContext& /*executionContext*/) const
{
  auto pOutputMode = filterArgs.value<ChoicesParameter::ValueType>(k_OutputMode_Key);
  auto pOutputUnits = filterArgs.value<ChoicesParameter::ValueType>(k_OutputUnits_Key);
  auto pApplySmoothing = filterArgs.value<bool>(k_ApplySmoothing_Key);
  auto pBinSizeDeg = filterArgs.value<float32>(k_BinSizeDeg_Key);
  auto pEulerAnglesPath = filterArgs.value<DataPath>(k_EulerAngles_Key);
  auto pPhasesPath = filterArgs.value<DataPath>(k_Phases_Key);
  auto pCrystalStructuresPath = filterArgs.value<DataPath>(k_CrystalStructures_Key);
  const auto materialNamesPath = filterArgs.value<DataPath>(k_MaterialNames_Key);
  const int32 selectedPhase = filterArgs.value<int32>(k_PhaseIndex_Key);
  auto pUseMask = filterArgs.value<bool>(k_UseMask_Key);
  auto pMaskPath = filterArgs.value<DataPath>(k_MaskArrayPath_Key);
  auto pOutputImageGeomPath = filterArgs.value<DataPath>(k_OutputImageGeometry_Key);
  auto pCellAttrMatName = filterArgs.value<DataObjectNameParameter::ValueType>(k_CellAttrMatName_Key);
  auto pExistingOdfGeomPath = filterArgs.value<DataPath>(k_ExistingOdfGeometry_Key);
  auto pComponentName = filterArgs.value<DataObjectNameParameter::ValueType>(k_ComponentName_Key);

  const bool isAppendMode = (pOutputMode == 1ULL);

  nx::core::Result<OutputActions> resultOutputActions;
  std::vector<PreflightValue> preflightUpdatedValues;

  const auto& crystalStructures = dataStructure.getDataRefAs<UInt32Array>(pCrystalStructuresPath);
  const auto* materialNamesPtr = dataStructure.getDataAs<StringArray>(materialNamesPath);
  if(materialNamesPtr == nullptr)
  {
    return {.outputActions =
                MakeErrorResult<OutputActions>(-12214, fmt::format("Material Names path '{}' must select a StringArray; the selected object has a different type.", materialNamesPath.toString()))};
  }
  const usize ensembleCount = crystalStructures.getNumberOfTuples();
  if(materialNamesPtr->getNumberOfTuples() != ensembleCount)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(
                -12215, fmt::format("Crystal Structures '{}' contains {} tuples, but Material Names '{}' contains {} tuples. Select arrays with equal ensemble tuple counts.",
                                    pCrystalStructuresPath.toString(), ensembleCount, materialNamesPath.toString(), materialNamesPtr->getNumberOfTuples()))};
  }
  if(selectedPhase <= 0 || static_cast<usize>(selectedPhase) >= ensembleCount)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12216, fmt::format("Phase Index ({}) must be positive and less than the ensemble tuple count ({}) in '{}' and '{}'.", selectedPhase,
                                                                                ensembleCount, pCrystalStructuresPath.toString(), materialNamesPath.toString()))};
  }
  // Preflight placeholders expose tuple counts but do not provide crystal values.
  if(crystalStructures.getIDataStoreRef().getStoreType() != IDataStore::StoreType::Empty)
  {
    const uint32 crystalCode = crystalStructures[static_cast<usize>(selectedPhase)];
    const auto allOps = ebsdlib::LaueOps::GetAllOrientationOps();
    if(crystalCode >= allOps.size() || allOps[crystalCode] == nullptr)
    {
      return {.outputActions =
                  MakeErrorResult<OutputActions>(-12217, fmt::format("Crystal Structures '{}' contains unsupported crystal-structure value {} at phase {}. Select a code from 0 through {}.",
                                                                     pCrystalStructuresPath.toString(), crystalCode, selectedPhase, allOps.size() - 1))};
    }
  }

  // Create mode uses the requested bin size. Append mode derives the bin size from the existing geometry's uniform spacing.
  double binSize = 0.0;
  int32_t nphi1 = 0;
  int32_t nPHI = 0;
  int32_t nphi2 = 0;
  bool binSizeIsDerived = false;

  if(!isAppendMode)
  {
    if(pBinSizeDeg <= 0.0f)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12200, fmt::format("Bin size must be > 0; got {}.", pBinSizeDeg))};
    }
    binSize = static_cast<double>(pBinSizeDeg);
    const double phi1Quotient = 360.0 / binSize;
    const double phiQuotient = 180.0 / binSize;
    constexpr double k_DivisorTolerance = 1.0e-6;
    if(std::abs(phi1Quotient - std::round(phi1Quotient)) > k_DivisorTolerance || std::abs(phiQuotient - std::round(phiQuotient)) > k_DivisorTolerance || std::round(phi1Quotient) <= 0.0 ||
       std::round(phiQuotient) <= 0.0)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12200, fmt::format("Bin size {} deg does not evenly divide 360 and 180. "
                                                                                  "Pick a value such that 360/binSize and 180/binSize are "
                                                                                  "both positive integers (e.g. 5.0, 2.0, 1.0).",
                                                                                  pBinSizeDeg))};
    }
    nphi1 = static_cast<int32_t>(std::round(360.0 / binSize));
    nPHI = static_cast<int32_t>(std::round(180.0 / binSize));
    nphi2 = nphi1;
  }
  else
  {
    // Append mode: validate the existing ODF ImageGeom up front so we can
    // derive the bin size.
    const auto* existingGeomPtr = dataStructure.getDataAs<ImageGeom>(pExistingOdfGeomPath);
    if(existingGeomPtr == nullptr)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12206, fmt::format("Existing ODF Image Geometry '{}' must exist and "
                                                                                  "be an ImageGeom.",
                                                                                  pExistingOdfGeomPath.toString()))};
    }
    const FloatVec3 spacing = existingGeomPtr->getSpacing();
    constexpr float32 k_SpacingTolerance = 1.0e-6f;
    if(std::abs(spacing[0] - spacing[1]) > k_SpacingTolerance || std::abs(spacing[0] - spacing[2]) > k_SpacingTolerance)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12207, fmt::format("Existing ODF Image Geometry '{}' has "
                                                                                  "non-uniform spacing ({}, {}, {}) across axes; "
                                                                                  "ODF requires uniform bins on all three axes.",
                                                                                  pExistingOdfGeomPath.toString(), spacing[0], spacing[1], spacing[2]))};
    }
    if(spacing[0] <= 0.0f)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12207, fmt::format("Existing ODF Image Geometry '{}' has zero or negative "
                                                                                  "spacing ({}); cannot derive a valid bin size.",
                                                                                  pExistingOdfGeomPath.toString(), spacing[0]))};
    }
    binSize = static_cast<double>(spacing[0]);
    binSizeIsDerived = true;
    // Derived dims from existing geometry: X = phi2, Y = PHI, Z = phi1.
    nphi2 = static_cast<int32_t>(existingGeomPtr->getNumXCells());
    nPHI = static_cast<int32_t>(existingGeomPtr->getNumYCells());
    nphi1 = static_cast<int32_t>(existingGeomPtr->getNumZCells());

    // Defense-in-depth: reject existing geometries whose dims do not match the
    // 360/binSize (phi1, phi2) and 180/binSize (PHI) integer-divisor
    // constraints. Catches both zero-sized axes and mismatched shape/spacing
    // combinations (e.g. dims 100x50x100 with spacing 5.0 deg).
    const double expected360 = 360.0 / binSize;
    const double expected180 = 180.0 / binSize;
    constexpr double k_DivisorTolerance = 1.0e-6;
    if(std::abs(expected360 - std::round(expected360)) > k_DivisorTolerance || std::abs(expected180 - std::round(expected180)) > k_DivisorTolerance ||
       std::round(expected360) != static_cast<double>(nphi1) || std::round(expected180) != static_cast<double>(nPHI) || std::round(expected360) != static_cast<double>(nphi2))
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12212, fmt::format("Existing ODF geometry is inconsistent with a "
                                                                                  "valid ODF grid: dims {}x{}x{} (phi1 x PHI x "
                                                                                  "phi2) with derived bin size {:.4f} deg "
                                                                                  "do not satisfy the 360/binSize (phi1/phi2) and "
                                                                                  "180/binSize (PHI) integer-divisor constraints.",
                                                                                  nphi1, nPHI, nphi2, binSize))};
    }
  }

  // Safety cap on total bin count. Each parallel worker allocates its own
  // std::vector<double> of this size, so a pathological user-supplied
  // bin_size_deg (e.g. 0.1 → 360x180x360 ~ 23M bins x 8 bytes x N threads) can
  // blow past usable RAM. 10^8 bins = 800 MB per accumulator is the hard upper
  // bound we allow.
  constexpr std::size_t k_MaxBinCount = 100000000;
  const std::size_t totalBinCount = static_cast<std::size_t>(nphi1) * static_cast<std::size_t>(nPHI) * static_cast<std::size_t>(nphi2);
  if(totalBinCount > k_MaxBinCount)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12211, fmt::format("Total bin count {} exceeds the safety limit of {}. Consider a "
                                                                                "larger bin_size_deg (e.g. 5.0 deg -> 72x36x72 = 186624 bins is "
                                                                                "typical; "
                                                                                "1.0 deg -> 360x180x360 ~ 23M bins is already above the cap).",
                                                                                totalBinCount, k_MaxBinCount))};
  }

  // Ensure all "Cell" arrays have the same number of tuples
  std::vector<DataPath> dataPaths;
  dataPaths.push_back(pPhasesPath);
  dataPaths.push_back(pEulerAnglesPath);
  if(pUseMask)
  {
    dataPaths.push_back(pMaskPath);
  }
  auto tupleValidityCheck = dataStructure.validateNumberOfTuples(dataPaths);
  if(!tupleValidityCheck)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-651, fmt::format("The following DataArrays all must have equal number "
                                                                              "of tuples but this was not satisfied.\n{}",
                                                                              tupleValidityCheck.error()))};
  }

  // 3. Same parent (AttributeMatrix) for the cell-level inputs.
  const DataPath eulerParent = pEulerAnglesPath.getParent();
  const DataPath phasesParent = pPhasesPath.getParent();
  if(eulerParent != phasesParent)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12205, fmt::format("Cell Euler Angles ('{}') and Cell Phases ('{}') must live "
                                                                                "on the same AttributeMatrix.",
                                                                                pEulerAnglesPath.toString(), pPhasesPath.toString()))};
  }
  if(pUseMask)
  {
    const DataPath maskParent = pMaskPath.getParent();
    if(maskParent != eulerParent)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12205, fmt::format("Mask ('{}') must live on the same AttributeMatrix as "
                                                                                  "the Cell Euler Angles ('{}').",
                                                                                  pMaskPath.toString(), pEulerAnglesPath.toString()))};
    }
  }

  // 4. Mode-specific: resolve target cell AttributeMatrix path + component
  // path, then check collisions.
  DataPath targetImageGeomPath;
  DataPath cellAttrMatPath;
  if(!isAppendMode)
  {
    targetImageGeomPath = pOutputImageGeomPath;
    cellAttrMatPath = pOutputImageGeomPath.createChildPath(pCellAttrMatName);
  }
  else
  {
    // Append mode: discover the existing geometry's cell AttributeMatrix via
    // getCellData()/getCellDataPath(). getCellDataPath() throws if there is no
    // cell-data AttributeMatrix assigned, so use the nullable accessor first
    // and fail gracefully with an error code instead.
    const auto& existingGeom = dataStructure.getDataRefAs<ImageGeom>(pExistingOdfGeomPath);
    if(existingGeom.getCellData() == nullptr)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12208, fmt::format("Existing ODF Image Geometry '{}' has no cell "
                                                                                  "AttributeMatrix; cannot append a component.",
                                                                                  pExistingOdfGeomPath.toString()))};
    }
    cellAttrMatPath = existingGeom.getCellDataPath();
    targetImageGeomPath = pExistingOdfGeomPath;

    // Component-name collision: the new component must not already exist on the
    // target cell AM.
    const DataPath proposedComponentPath = cellAttrMatPath.createChildPath(pComponentName);
    if(dataStructure.getDataAs<IDataArray>(proposedComponentPath) != nullptr)
    {
      return {.outputActions = MakeErrorResult<OutputActions>(-12209, fmt::format("Component '{}' already exists on the target geometry's cell "
                                                                                  "AttributeMatrix ('{}'); pick a different component name.",
                                                                                  pComponentName, cellAttrMatPath.toString()))};
    }
  }

  const DataPath metadataGroupPath = targetImageGeomPath.createChildPath(mtrsim::odf_metadata::k_GroupName.str());
  const DataPath unitsPath = mtrsim::odf_metadata::CreateUnitsPath(targetImageGeomPath, pComponentName);
  const auto* metadataObjectPtr = dataStructure.getData(metadataGroupPath);
  if((!isAppendMode && cellAttrMatPath == metadataGroupPath) || (metadataObjectPtr != nullptr && dataStructure.getDataAs<DataGroup>(metadataGroupPath) == nullptr))
  {
    return {.outputActions = MakeErrorResult<OutputActions>(
                -12218, fmt::format("ODF metadata path '{}' must be a DataGroup. Choose a different cell Attribute Matrix name or remove the conflicting object.", metadataGroupPath.toString()))};
  }
  if(dataStructure.getData(unitsPath) != nullptr)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12219, fmt::format("ODF units metadata '{}' already exists. Choose a different component name.", unitsPath.toString()))};
  }

  // 5. Output actions.
  if(!isAppendMode)
  {
    // ImageGeom XYZ = (phi2, PHI, phi1); array tuple shape ZYX = (nphi1, nPHI,
    // nphi2).
    const std::vector<usize> imageGeomDimsXYZ = {static_cast<usize>(nphi2), static_cast<usize>(nPHI), static_cast<usize>(nphi1)};
    const std::vector<float32> imageGeomOrigin = {0.0f, 0.0f, 0.0f};
    const std::vector<float32> imageGeomSpacingXYZ = {static_cast<float32>(binSize), static_cast<float32>(binSize), static_cast<float32>(binSize)};
    auto createImageGeomAction = std::make_unique<CreateImageGeometryAction>(pOutputImageGeomPath, imageGeomDimsXYZ, imageGeomOrigin, imageGeomSpacingXYZ, pCellAttrMatName);
    resultOutputActions.value().appendAction(std::move(createImageGeomAction));
  }
  {
    const std::vector<usize> arrayTupleShapeZYX = {static_cast<usize>(nphi1), static_cast<usize>(nPHI), static_cast<usize>(nphi2)};
    const std::vector<usize> componentShape = {1};
    const DataPath componentPath = cellAttrMatPath.createChildPath(pComponentName);
    auto createArrayAction = std::make_unique<CreateArrayAction>(DataType::float64, arrayTupleShapeZYX, componentShape, componentPath);
    resultOutputActions.value().appendAction(std::move(createArrayAction));
  }

  if(metadataObjectPtr == nullptr)
  {
    resultOutputActions.value().appendAction(std::make_unique<CreateDataGroupAction>(metadataGroupPath));
  }
  resultOutputActions.value().appendAction(std::make_unique<CreateStringArrayAction>(std::vector<usize>{1}, unitsPath, mtrsim::odf_metadata::UnitsToString(pOutputUnits)));

  // 6. PreflightUpdatedValues so the UI can preview the dims/settings.
  const std::string binSizeLabel = binSizeIsDerived ? fmt::format("{:.4f} (derived from existing geometry; user input ignored)", binSize) : fmt::format("{:.4f}", binSize);
  preflightUpdatedValues.push_back({"Output Mode", isAppendMode ? std::string{"Append to Existing ODF Geometry"} : std::string{"Create New ODF Geometry"}});
  preflightUpdatedValues.push_back({"Target ODF Image Geometry", targetImageGeomPath.toString()});
  preflightUpdatedValues.push_back({"Bin Size [deg]", binSizeLabel});
  preflightUpdatedValues.push_back({"Bins (phi1 x PHI x phi2)", fmt::format("{} x {} x {}", nphi1, nPHI, nphi2)});
  preflightUpdatedValues.push_back({"ImageGeom Dimensions (X, Y, Z)", fmt::format("{} x {} x {}", nphi2, nPHI, nphi1)});
  preflightUpdatedValues.push_back({"Total Bin Count", std::to_string(totalBinCount)});
  preflightUpdatedValues.push_back({"Smoothing", pApplySmoothing ? std::string{"enabled"} : std::string{"disabled"}});
  preflightUpdatedValues.push_back({"Output Units", (pOutputUnits == 1ULL) ? std::string{"MUD (Jacobian sin(PHI) applied)"} : std::string{"Count-Density (raw histogram)"}});
  preflightUpdatedValues.push_back({"Mask", pUseMask ? pMaskPath.toString() : std::string{"(none)"}});
  preflightUpdatedValues.push_back({"Selected Phase", std::to_string(selectedPhase)});

  const bool writePng = filterArgs.value<bool>(k_WritePng_Key);
  const bool saveAsImageGeometry = filterArgs.value<bool>(k_SaveAsImageGeometry_Key);
  if(!writePng && !saveAsImageGeometry)
  {
    return {.outputActions = std::move(resultOutputActions), .outputValues = std::move(preflightUpdatedValues)};
  }
  const int32 sectionCount = filterArgs.value<int32>(k_NumSections_Key);
  const int32 sectionsPerRow = filterArgs.value<int32>(k_SectionsPerRow_Key);
  const int32 sectionWidth = filterArgs.value<int32>(k_SectionWidth_Key);
  const auto colorTablePreset = filterArgs.value<CreateColorMapParameter::ValueType>(k_ColorTablePreset_Key);
  const auto scaleMode = filterArgs.value<ChoicesParameter::ValueType>(k_ScaleMode_Key);
  const float64 manualMinimum = filterArgs.value<float64>(k_ManualMinimum_Key);
  const float64 manualMaximum = filterArgs.value<float64>(k_ManualMaximum_Key);
  const auto outputPngPath = filterArgs.value<FileSystemPathParameter::ValueType>(k_OutputPngPath_Key);
  const auto outputImageGeometryPath = filterArgs.value<DataPath>(k_OutputPresentationGeometry_Key);

  const bool deferredMetadata = crystalStructures.getIDataStoreRef().getStoreType() == IDataStore::StoreType::Empty;
  const bool hasMaterialName = !deferredMetadata && !materialNamesPtr->isPlaceholder();
  // Empty ensemble stores provide no values. Triclinic limits reserve the largest presentation page.
  const uint32 crystalStructure = deferredMetadata ? ebsdlib::CrystalStructure::Triclinic : crystalStructures[static_cast<usize>(selectedPhase)];
  const auto limits = ebsdlib::GetODFEulerPlotLimits(crystalStructure);
  if(sectionCount < 2)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12309, fmt::format("Number of Sections ({}) must be at least 2.", sectionCount))};
  }
  if(sectionsPerRow < 1 || sectionsPerRow > sectionCount)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12310, fmt::format("Sections Per Row ({}) must be in [1, Number of Sections ({})].", sectionsPerRow, sectionCount))};
  }
  if(sectionWidth < 128)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12311, fmt::format("Section Width ({}) pixels must be at least 128 pixels.", sectionWidth))};
  }

  auto colorControls = ColorTableUtilities::ExtractControlPoints(colorTablePreset);
  if(colorControls.invalid())
  {
    return {.outputActions = MakeErrorResult<OutputActions>(-12312, fmt::format("Color Table Preset '{}' is not a valid RGB preset.", colorTablePreset))};
  }
  if(scaleMode == 1 && (!std::isfinite(manualMinimum) || !std::isfinite(manualMaximum) || manualMinimum < 0.0 || manualMinimum >= manualMaximum))
  {
    return {.outputActions = MakeErrorResult<OutputActions>(
                -12313, fmt::format("Manual MUD range must have finite endpoints and satisfy 0 <= minimum < maximum; received [{}, {}].", manualMinimum, manualMaximum))};
  }
  if(writePng && outputPngPath.extension() != ".png")
  {
    return {.outputActions =
                MakeErrorResult<OutputActions>(-12314, fmt::format("Output PNG file '{}' has extension '{}'; use the '.png' extension.", outputPngPath.string(), outputPngPath.extension().string()))};
  }

  ebsdlib::ODFSectionConfiguration config;
  // Image geometry XYZ corresponds to phi2, PHI, phi1; the borrowed values retain phi2-fastest storage.
  config.grid.dimensions = {static_cast<usize>(nphi1), static_cast<usize>(nPHI), static_cast<usize>(nphi2)};
  config.grid.originDeg = {0.0, 0.0, 0.0};
  config.grid.spacingDeg = {binSize, binSize, binSize};
  config.grid.laueOpsIndex = crystalStructure;
  config.sectionCount = static_cast<usize>(sectionCount);
  config.sectionsPerRow = static_cast<usize>(sectionsPerRow);
  config.sectionWidth = sectionWidth;
  config.scaleMode = scaleMode == 0 ? ebsdlib::ODFScaleMode::Automatic : ebsdlib::ODFScaleMode::Manual;
  config.manualMinimumMUD = manualMinimum;
  config.manualMaximumMUD = manualMaximum;
  config.colorControlPoints = std::move(colorControls.value());
  config.title = filterArgs.value<StringParameter::ValueType>(k_ODF_Image_Title_Key);
  if(hasMaterialName)
  {
    config.materialName = (*materialNamesPtr)[static_cast<usize>(selectedPhase)];
  }
  config.phaseNumber = selectedPhase;

  const auto layout = ebsdlib::ComputeODFSectionLayout(config);
  if(layout.pageWidth <= 0 || layout.pageHeight <= 0)
  {
    return {.outputActions = MakeErrorResult<OutputActions>(
                -12315, fmt::format("ODF section layout is not representable for Section Width {}, Number of Sections {}, and Sections Per Row {}.", sectionWidth, sectionCount, sectionsPerRow))};
  }

  if(saveAsImageGeometry)
  {
    resultOutputActions.value().appendAction(
        std::make_unique<CreateImageGeometryAction>(outputImageGeometryPath, CreateImageGeometryAction::DimensionType{static_cast<usize>(layout.pageWidth), static_cast<usize>(layout.pageHeight), 1},
                                                    CreateImageGeometryAction::OriginType{0.0F, 0.0F, 0.0F}, CreateImageGeometryAction::SpacingType{1.0F, 1.0F, 1.0F}, std::string{"Cell Data"}));
    const DataPath rgbArrayPath = outputImageGeometryPath.createChildPath(std::string{"Cell Data"}).createChildPath(std::string{"ODF Sections RGB"});
    resultOutputActions.value().appendAction(
        std::make_unique<CreateArrayAction>(DataType::uint8, ShapeType{1, static_cast<usize>(layout.pageHeight), static_cast<usize>(layout.pageWidth)}, ShapeType{3}, rgbArrayPath));
  }

  if(!deferredMetadata)
  {
    const auto sectionAngles = ebsdlib::GenerateODFSectionAngles(crystalStructure, config.sectionCount);
    preflightUpdatedValues.push_back({"Section Angles", fmt::format("{} degrees", fmt::join(sectionAngles, ", "))});
    if(hasMaterialName)
    {
      preflightUpdatedValues.push_back({"Material", config.materialName});
    }
    preflightUpdatedValues.push_back({"Laue Class", ebsdlib::LaueOps::GetLaueNames()[crystalStructure]});
    preflightUpdatedValues.push_back({"Euler Plot Limits", fmt::format("phi1 < {}, PHI < {}, phi2 < {} degrees", limits.phi1MaxDeg, limits.phiMaxDeg, limits.phi2MaxDeg)});
  }
  preflightUpdatedValues.push_back({"Input Units", mtrsim::odf_metadata::UnitsToString(pOutputUnits)});
  preflightUpdatedValues.push_back({"Rendered Units", "MUD"});
  preflightUpdatedValues.push_back({"Full Grid Dimensions (phi1 x PHI x phi2)", fmt::format("{} x {} x {}", nphi1, nPHI, nphi2)});
  preflightUpdatedValues.push_back({"Full Grid Spacing (phi1, PHI, phi2)", fmt::format("{0}, {0}, {0} degrees", binSize)});
  preflightUpdatedValues.push_back({"Section Layout", fmt::format("{} columns x {} rows", layout.columns, layout.rows)});
  preflightUpdatedValues.push_back(
      {"Composite Dimensions", deferredMetadata ? fmt::format("{} x {} pixels (conservative maximum estimate; exact dimensions are set during execution).", layout.pageWidth, layout.pageHeight) :
                                                  fmt::format("{} x {} pixels", layout.pageWidth, layout.pageHeight)});
  preflightUpdatedValues.push_back({"Scale Mode", scaleMode == 0 ? std::string{"Automatic"} : fmt::format("Manual: [{}, {}] MUD", manualMinimum, manualMaximum)});

  return {.outputActions = std::move(resultOutputActions), .outputValues = std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> ComputeODFFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* /*pipelineNode*/, const MessageHandler& messageHandler,
                                       const std::atomic_bool& shouldCancel, const ExecutionContext& /*executionContext*/) const
{
  const auto pOutputMode = filterArgs.value<ChoicesParameter::ValueType>(k_OutputMode_Key);
  const bool isAppendMode = (pOutputMode == 1ULL);

  ComputeODFInputValues inputValues;
  usize contributingVoxelCount = 0;
  inputValues.contributingVoxelCount = &contributingVoxelCount;

  inputValues.selectedPhase = filterArgs.value<int32>(k_PhaseIndex_Key);
  inputValues.applySmoothing = filterArgs.value<bool>(k_ApplySmoothing_Key);
  inputValues.outputUnits = filterArgs.value<ChoicesParameter::ValueType>(k_OutputUnits_Key);
  inputValues.useMask = filterArgs.value<bool>(k_UseMask_Key);
  inputValues.eulerAnglesPath = filterArgs.value<DataPath>(k_EulerAngles_Key);
  inputValues.phasesPath = filterArgs.value<DataPath>(k_Phases_Key);
  inputValues.crystalStructuresPath = filterArgs.value<DataPath>(k_CrystalStructures_Key);
  inputValues.maskPath = filterArgs.value<DataPath>(k_MaskArrayPath_Key);
  inputValues.componentName = filterArgs.value<DataObjectNameParameter::ValueType>(k_ComponentName_Key);

  if(!isAppendMode)
  {
    inputValues.binSizeDeg = static_cast<double>(filterArgs.value<float32>(k_BinSizeDeg_Key));
    inputValues.outputImageGeometry = filterArgs.value<DataPath>(k_OutputImageGeometry_Key);
    inputValues.cellAttrMatName = filterArgs.value<DataObjectNameParameter::ValueType>(k_CellAttrMatName_Key);

    // Pre-compute axis dims here so the algorithm does not have to repeat it.
    // Preflight has already validated that binSizeDeg cleanly divides 360 and
    // 180, so std::round is exact.
    inputValues.nphi1 = static_cast<int32_t>(std::round(360.0 / inputValues.binSizeDeg));
    inputValues.nPHI = static_cast<int32_t>(std::round(180.0 / inputValues.binSizeDeg));
    inputValues.nphi2 = inputValues.nphi1;
  }
  else
  {
    // Append mode: override target geometry / cell attr matrix name / bin size
    // from the existing geometry.
    const auto pExistingOdfGeomPath = filterArgs.value<DataPath>(k_ExistingOdfGeometry_Key);
    const auto& existingGeom = dataStructure.getDataRefAs<ImageGeom>(pExistingOdfGeomPath);
    const FloatVec3 spacing = existingGeom.getSpacing();

    inputValues.outputImageGeometry = pExistingOdfGeomPath;
    inputValues.cellAttrMatName = existingGeom.getCellDataPath().getTargetName();
    inputValues.binSizeDeg = static_cast<double>(spacing[0]);
    inputValues.nphi2 = static_cast<int32_t>(existingGeom.getNumXCells());
    inputValues.nPHI = static_cast<int32_t>(existingGeom.getNumYCells());
    inputValues.nphi1 = static_cast<int32_t>(existingGeom.getNumZCells());
  }

  auto result = ComputeODF(dataStructure, messageHandler, shouldCancel, &inputValues)();
  if(result.invalid() || shouldCancel || (!filterArgs.value<bool>(k_WritePng_Key) && !filterArgs.value<bool>(k_SaveAsImageGeometry_Key)))
  {
    return result;
  }

  if(contributingVoxelCount == 0)
  {
    return MergeResults(
        std::move(result),
        MakeErrorResult(-12324, fmt::format("Cannot create ODF presentation for phase {}: no voxels contributed. Check the phase labels and mask, or disable both presentation outputs.",
                                            inputValues.selectedPhase)));
  }

  const auto materialNamesPath = filterArgs.value<DataPath>(k_MaterialNames_Key);
  const auto& materialNames = dataStructure.getDataRefAs<StringArray>(materialNamesPath);
  if(materialNames.isPlaceholder())
  {
    return MergeResults(std::move(result),
                        MakeErrorResult(-12323, fmt::format("Material Names '{}' values for phase {} are not available during ODF presentation. Execute the upstream material-data reader first.",
                                                            materialNamesPath.toString(), inputValues.selectedPhase)));
  }

  GenerateODFSectionsInputValues presentationValues;
  presentationValues.inputImageGeometryPath = inputValues.outputImageGeometry;
  presentationValues.odfArrayPath = inputValues.outputImageGeometry.createChildPath(inputValues.cellAttrMatName).createChildPath(inputValues.componentName);
  presentationValues.units = inputValues.outputUnits == 0 ? ebsdlib::ODFValueUnits::CountDensity : ebsdlib::ODFValueUnits::MUD;
  presentationValues.phaseIndex = inputValues.selectedPhase;
  presentationValues.crystalStructure = dataStructure.getDataRefAs<UInt32Array>(inputValues.crystalStructuresPath)[static_cast<usize>(inputValues.selectedPhase)];
  presentationValues.materialName = materialNames[static_cast<usize>(inputValues.selectedPhase)];
  presentationValues.title = filterArgs.value<StringParameter::ValueType>(k_ODF_Image_Title_Key);
  presentationValues.numberOfSections = filterArgs.value<int32>(k_NumSections_Key);
  presentationValues.sectionsPerRow = filterArgs.value<int32>(k_SectionsPerRow_Key);
  presentationValues.sectionWidth = filterArgs.value<int32>(k_SectionWidth_Key);
  presentationValues.colorTablePreset = filterArgs.value<CreateColorMapParameter::ValueType>(k_ColorTablePreset_Key);
  presentationValues.scaleMode = filterArgs.value<ChoicesParameter::ValueType>(k_ScaleMode_Key) == 0 ? ebsdlib::ODFScaleMode::Automatic : ebsdlib::ODFScaleMode::Manual;
  presentationValues.manualMinimumMUD = filterArgs.value<float64>(k_ManualMinimum_Key);
  presentationValues.manualMaximumMUD = filterArgs.value<float64>(k_ManualMaximum_Key);
  presentationValues.writePng = filterArgs.value<bool>(k_WritePng_Key);
  presentationValues.outputPngPath = filterArgs.value<FileSystemPathParameter::ValueType>(k_OutputPngPath_Key);
  presentationValues.saveAsImageGeometry = filterArgs.value<bool>(k_SaveAsImageGeometry_Key);
  presentationValues.outputImageGeometryPath = filterArgs.value<DataPath>(k_OutputPresentationGeometry_Key);

  return MergeResults(std::move(result), GenerateODFSections(dataStructure, messageHandler, shouldCancel, &presentationValues)());
}

//------------------------------------------------------------------------------
Result<Arguments> ComputeODFFilter::FromSIMPLJson(const nlohmann::json& /*json*/)
{
  Arguments args = ComputeODFFilter().getDefaultArguments();

  std::vector<Result<>> results;

  /* This is a NEW filter and has no SIMPL (DREAM3D v6) equivalent. */

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult), std::move(args));
}

} // namespace nx::core

#include "ComputeODFFilter.hpp"

#include "MTRSim/Filters/Algorithms/ComputeODF.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/Geometry/IGeometry.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/IDataArray.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"

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
  params.insert(std::make_unique<BoolParameter>(k_ApplySmoothing_Key, "Apply Smoothing",
                                                "When true, each symmetric Euler tuple is distributed across 27 bins (1 center + 6 faces + 12 edges + 8 corners) using "
                                                "the MATLAB calc_ODF.m tri-linear weights. When false, each tuple deposits 1.0 into its center bin only.",
                                                true));
  params.insert(std::make_unique<Float32Parameter>(k_BinSizeDeg_Key, "Bin Size [deg]",
                                                   "Uniform bin size in degrees, used for all three Bunge-Euler axes (phi1, PHI, phi2). 360/binSize and 180/binSize "
                                                   "must both yield a positive integer (e.g. 5.0 -> 72 x 36 x 72 bins).",
                                                   5.0f));
  params.insert(std::make_unique<ChoicesParameter>(k_OutputUnits_Key, "Output Units",
                                                   "Count-Density emits the raw normalized histogram (sums to 1 over all Bunge bins; bit-exact match to MATLAB calc_ODF.m). "
                                                   "MUD applies the per-PHI-row sin(PHI) Jacobian correction so values represent density on SO(3) divided by uniform "
                                                   "density (matches MTEX plot(odf) convention; sub-uniform regions read < 1, peak texture reads several MUD).",
                                                   0ULL, ChoicesParameter::Choices{"Count-Density", "MUD"}));

  params.insertSeparator(Parameters::Separator{"Required Input Cell Data"});
  params.insert(std::make_unique<ArraySelectionParameter>(k_EulerAngles_Key, "Cell Euler Angles",
                                                          "Per-voxel Bunge Euler angles (phi1, PHI, phi2) in radians. Must be a 3-component Float32 array on a cell-level "
                                                          "AttributeMatrix.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::float32}, ArraySelectionParameter::AllowedComponentShapes{{3}}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_Phases_Key, "Cell Phases",
                                                          "Per-voxel phase label (1-based; 0 marks excluded voxels). Must be a 1-component Int32 array on the same "
                                                          "AttributeMatrix as the Euler angles.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::int32}, ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insertSeparator(Parameters::Separator{"Optional Input Cell Data"});
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_UseMask_Key, "Use Mask Array",
                                                                 "When enabled, only voxels where the selected mask is true (and phase != 0) contribute to the ODF.", false));
  params.insert(std::make_unique<ArraySelectionParameter>(k_Mask_Key, "Mask",
                                                          "Per-voxel boolean mask. Must be a 1-component Bool array on the same AttributeMatrix as the Euler angles. Only "
                                                          "consulted when 'Use Mask Array' is enabled.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::boolean, DataType::uint8}, ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insertSeparator(Parameters::Separator{"Required Input Ensemble Data"});
  params.insert(std::make_unique<ArraySelectionParameter>(k_CrystalStructures_Key, "Crystal Structures",
                                                          "Per-phase EbsdLib crystal-structure code (e.g. 0 = Hexagonal_High, 1 = Cubic_High). Must be a 1-component "
                                                          "UInt32 array on the ensemble-level AttributeMatrix; indexed by phase label.",
                                                          DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::uint32}, ArraySelectionParameter::AllowedComponentShapes{{1}}));


  params.insertSeparator(Parameters::Separator{"Output Mode"});
  params.insertLinkableParameter(std::make_unique<ChoicesParameter>(k_OutputMode_Key, "Output Mode",
                                                                    "Whether to create a new ODF ImageGeometry or append a new component to an existing one. In Append mode the bin "
                                                                    "size is derived from the existing geometry's spacing and the user's bin size input is ignored.",
                                                                    0ULL, ChoicesParameter::Choices{"Create New ODF Geometry", "Append to Existing ODF Geometry"}));

  params.insertSeparator(Parameters::Separator{"Output Data Object(s)"});
  params.insert(std::make_unique<DataGroupCreationParameter>(k_OutputImageGeometry_Key, "Output ODF Image Geometry",
                                                             "Path at which the ODF ImageGeom will be created. Its (X, Y, Z) dimensions map to (phi2, PHI, phi1).", DataPath({"ODF"})));
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellAttrMatName_Key, "Cell Attribute Matrix Name",
                                                          "Name of the cell AttributeMatrix created under the output ImageGeom. The Float64 ODF array is placed inside it.", "Cell Data"));
  params.insert(std::make_unique<GeometrySelectionParameter>(k_ExistingOdfGeometry_Key, "Existing ODF Image Geometry",
                                                             "ImageGeom representing an existing Euler-space ODF grid that this filter will append a new component to. Spacing must be "
                                                             "uniform across all three axes.",
                                                             DataPath{}, GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<DataObjectNameParameter>(k_ComponentName_Key, "ODF Component Array Name",
                                                          "Name of the Float64 single-component ODF array created on the output cell AttributeMatrix.", "Component 1"));

  params.linkParameters(k_UseMask_Key, k_Mask_Key, true);
  params.linkParameters(k_OutputMode_Key, k_OutputImageGeometry_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  params.linkParameters(k_OutputMode_Key, k_CellAttrMatName_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  params.linkParameters(k_OutputMode_Key, k_BinSizeDeg_Key, std::make_any<ChoicesParameter::ValueType>(0ULL));
  params.linkParameters(k_OutputMode_Key, k_ExistingOdfGeometry_Key, std::make_any<ChoicesParameter::ValueType>(1ULL));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ComputeODFFilter::parametersVersion() const
{
  return 3;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ComputeODFFilter::clone() const
{
  return std::make_unique<ComputeODFFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ComputeODFFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                         const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  auto pOutputMode = filterArgs.value<ChoicesParameter::ValueType>(k_OutputMode_Key);
  auto pOutputUnits = filterArgs.value<ChoicesParameter::ValueType>(k_OutputUnits_Key);
  auto pApplySmoothing = filterArgs.value<bool>(k_ApplySmoothing_Key);
  auto pBinSizeDeg = filterArgs.value<float32>(k_BinSizeDeg_Key);
  auto pEulerAnglesPath = filterArgs.value<DataPath>(k_EulerAngles_Key);
  auto pPhasesPath = filterArgs.value<DataPath>(k_Phases_Key);
  auto pCrystalStructuresPath = filterArgs.value<DataPath>(k_CrystalStructures_Key);
  auto pUseMask = filterArgs.value<bool>(k_UseMask_Key);
  auto pMaskPath = filterArgs.value<DataPath>(k_Mask_Key);
  auto pOutputImageGeomPath = filterArgs.value<DataPath>(k_OutputImageGeometry_Key);
  auto pCellAttrMatName = filterArgs.value<DataObjectNameParameter::ValueType>(k_CellAttrMatName_Key);
  auto pExistingOdfGeomPath = filterArgs.value<DataPath>(k_ExistingOdfGeometry_Key);
  auto pComponentName = filterArgs.value<DataObjectNameParameter::ValueType>(k_ComponentName_Key);

  const bool isAppendMode = (pOutputMode == 1ULL);

  nx::core::Result<OutputActions> resultOutputActions;
  std::vector<PreflightValue> preflightUpdatedValues;

  // 1. Determine bin size (& derived dims). In Create New mode we validate the user's input;
  //    in Append mode we derive from the existing geometry's (uniform) spacing and silently
  //    ignore the user's bin_size_deg input.
  double binSize = 0.0;
  int32_t nphi1 = 0;
  int32_t nPHI = 0;
  int32_t nphi2 = 0;
  bool binSizeIsDerived = false;

  if(!isAppendMode)
  {
    if(pBinSizeDeg <= 0.0f)
    {
      return {MakeErrorResult<OutputActions>(-12200, fmt::format("Bin size must be > 0; got {}.", pBinSizeDeg))};
    }
    binSize = static_cast<double>(pBinSizeDeg);
    const double phi1Quotient = 360.0 / binSize;
    const double phiQuotient = 180.0 / binSize;
    constexpr double k_DivisorTolerance = 1.0e-6;
    if(std::abs(phi1Quotient - std::round(phi1Quotient)) > k_DivisorTolerance || std::abs(phiQuotient - std::round(phiQuotient)) > k_DivisorTolerance || std::round(phi1Quotient) <= 0.0
       || std::round(phiQuotient) <= 0.0)
    {
      return {MakeErrorResult<OutputActions>(
          -12200,
          fmt::format("Bin size {} deg does not evenly divide 360 and 180. Pick a value such that 360/binSize and 180/binSize are both positive integers (e.g. 5.0, 2.0, 1.0).", pBinSizeDeg))};
    }
    nphi1 = static_cast<int32_t>(std::round(360.0 / binSize));
    nPHI = static_cast<int32_t>(std::round(180.0 / binSize));
    nphi2 = nphi1;
  }
  else
  {
    // Append mode: validate the existing ODF ImageGeom up front so we can derive the bin size.
    const auto* existingGeom = dataStructure.getDataAs<ImageGeom>(pExistingOdfGeomPath);
    if(existingGeom == nullptr)
    {
      return {MakeErrorResult<OutputActions>(-12206, fmt::format("Existing ODF Image Geometry '{}' must exist and be an ImageGeom.", pExistingOdfGeomPath.toString()))};
    }
    const FloatVec3 spacing = existingGeom->getSpacing();
    constexpr float32 k_SpacingTolerance = 1.0e-6f;
    if(std::abs(spacing[0] - spacing[1]) > k_SpacingTolerance || std::abs(spacing[0] - spacing[2]) > k_SpacingTolerance)
    {
      return {MakeErrorResult<OutputActions>(-12207,
                                             fmt::format("Existing ODF Image Geometry '{}' has non-uniform spacing ({}, {}, {}) across axes; ODF requires uniform bins on all three axes.",
                                                         pExistingOdfGeomPath.toString(), spacing[0], spacing[1], spacing[2]))};
    }
    if(spacing[0] <= 0.0f)
    {
      return {MakeErrorResult<OutputActions>(
          -12207, fmt::format("Existing ODF Image Geometry '{}' has zero or negative spacing ({}); cannot derive a valid bin size.", pExistingOdfGeomPath.toString(), spacing[0]))};
    }
    binSize = static_cast<double>(spacing[0]);
    binSizeIsDerived = true;
    // Derived dims from existing geometry: X = phi2, Y = PHI, Z = phi1.
    nphi2 = static_cast<int32_t>(existingGeom->getNumXCells());
    nPHI = static_cast<int32_t>(existingGeom->getNumYCells());
    nphi1 = static_cast<int32_t>(existingGeom->getNumZCells());

    // Defense-in-depth: reject existing geometries whose dims don't match the 360/binSize (phi1, phi2)
    // and 180/binSize (PHI) integer-divisor constraints. Catches both zero-sized axes and mismatched
    // shape/spacing combinations (e.g. dims 100x50x100 with spacing 5.0 deg).
    const double expected360 = 360.0 / binSize;
    const double expected180 = 180.0 / binSize;
    constexpr double k_DivisorTolerance = 1.0e-6;
    if(std::abs(expected360 - std::round(expected360)) > k_DivisorTolerance || std::abs(expected180 - std::round(expected180)) > k_DivisorTolerance
       || std::round(expected360) != static_cast<double>(nphi1) || std::round(expected180) != static_cast<double>(nPHI) || std::round(expected360) != static_cast<double>(nphi2))
    {
      return {MakeErrorResult<OutputActions>(-12212,
                                             fmt::format("Existing ODF geometry is inconsistent with a valid ODF grid: dims {}x{}x{} (phi1 x PHI x phi2) with derived bin size {:.4f} deg "
                                                         "do not satisfy the 360/binSize (phi1/phi2) and 180/binSize (PHI) integer-divisor constraints.",
                                                         nphi1, nPHI, nphi2, binSize))};
    }
  }

  // Safety cap on total bin count. Each parallel worker allocates its own std::vector<double>
  // of this size, so a pathological user-supplied bin_size_deg (e.g. 0.1 → 360x180x360 ~ 23M bins
  // x 8 bytes x N threads) can blow past usable RAM. 10^8 bins = 800 MB per accumulator is the
  // hard upper bound we allow.
  constexpr std::size_t k_MaxBinCount = 100000000;
  const std::size_t totalBinCount = static_cast<std::size_t>(nphi1) * static_cast<std::size_t>(nPHI) * static_cast<std::size_t>(nphi2);
  if(totalBinCount > k_MaxBinCount)
  {
    return {MakeErrorResult<OutputActions>(-12211,
                                           fmt::format("Total bin count {} exceeds the safety limit of {}. Consider a larger bin_size_deg (e.g. 5.0 deg -> 72x36x72 = 186624 bins is typical; "
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
    return {MakeErrorResult<OutputActions>(-651, fmt::format("The following DataArrays all must have equal number of tuples but this was not satisfied.\n{}", tupleValidityCheck.error()))};
  }


  // 3. Same parent (AttributeMatrix) for the cell-level inputs.
  const DataPath eulerParent = pEulerAnglesPath.getParent();
  const DataPath phasesParent = pPhasesPath.getParent();
  if(eulerParent != phasesParent)
  {
    return {MakeErrorResult<OutputActions>(
        -12205, fmt::format("Cell Euler Angles ('{}') and Cell Phases ('{}') must live on the same AttributeMatrix.", pEulerAnglesPath.toString(), pPhasesPath.toString()))};
  }
  if(pUseMask)
  {
    const DataPath maskParent = pMaskPath.getParent();
    if(maskParent != eulerParent)
    {
      return {
          MakeErrorResult<OutputActions>(-12205, fmt::format("Mask ('{}') must live on the same AttributeMatrix as the Cell Euler Angles ('{}').", pMaskPath.toString(), pEulerAnglesPath.toString()))};
    }
  }

  // 4. Mode-specific: resolve target cell AttributeMatrix path + component path, then check collisions.
  DataPath targetImageGeomPath;
  DataPath cellAttrMatPath;
  if(!isAppendMode)
  {
    targetImageGeomPath = pOutputImageGeomPath;
    cellAttrMatPath = pOutputImageGeomPath.createChildPath(pCellAttrMatName);
  }
  else
  {
    // Append mode: discover the existing geometry's cell AttributeMatrix via getCellData()/getCellDataPath().
    // getCellDataPath() throws if there's no cell-data AttributeMatrix assigned, so use the nullable
    // accessor first and fail gracefully with an error code instead.
    const auto& existingGeom = dataStructure.getDataRefAs<ImageGeom>(pExistingOdfGeomPath);
    if(existingGeom.getCellData() == nullptr)
    {
      return {MakeErrorResult<OutputActions>(
          -12208, fmt::format("Existing ODF Image Geometry '{}' has no cell AttributeMatrix; cannot append a component.", pExistingOdfGeomPath.toString()))};
    }
    cellAttrMatPath = existingGeom.getCellDataPath();
    targetImageGeomPath = pExistingOdfGeomPath;

    // Component-name collision: the new component must not already exist on the target cell AM.
    const DataPath proposedComponentPath = cellAttrMatPath.createChildPath(pComponentName);
    if(dataStructure.getDataAs<IDataArray>(proposedComponentPath) != nullptr)
    {
      return {MakeErrorResult<OutputActions>(
          -12209, fmt::format("Component '{}' already exists on the target geometry's cell AttributeMatrix ('{}'); pick a different component name.", pComponentName, cellAttrMatPath.toString()))};
    }
  }

  // 5. Output actions.
  if(!isAppendMode)
  {
    // ImageGeom XYZ = (phi2, PHI, phi1); array tuple shape ZYX = (nphi1, nPHI, nphi2).
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

  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> ComputeODFFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                       const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  const auto pOutputMode = filterArgs.value<ChoicesParameter::ValueType>(k_OutputMode_Key);
  const bool isAppendMode = (pOutputMode == 1ULL);

  ComputeODFInputValues inputValues;

  inputValues.applySmoothing = filterArgs.value<bool>(k_ApplySmoothing_Key);
  inputValues.outputUnits = filterArgs.value<ChoicesParameter::ValueType>(k_OutputUnits_Key);
  inputValues.useMask = filterArgs.value<bool>(k_UseMask_Key);
  inputValues.eulerAnglesPath = filterArgs.value<DataPath>(k_EulerAngles_Key);
  inputValues.phasesPath = filterArgs.value<DataPath>(k_Phases_Key);
  inputValues.crystalStructuresPath = filterArgs.value<DataPath>(k_CrystalStructures_Key);
  inputValues.maskPath = filterArgs.value<DataPath>(k_Mask_Key);
  inputValues.componentName = filterArgs.value<DataObjectNameParameter::ValueType>(k_ComponentName_Key);

  if(!isAppendMode)
  {
    inputValues.binSizeDeg = static_cast<double>(filterArgs.value<float32>(k_BinSizeDeg_Key));
    inputValues.outputImageGeometry = filterArgs.value<DataPath>(k_OutputImageGeometry_Key);
    inputValues.cellAttrMatName = filterArgs.value<DataObjectNameParameter::ValueType>(k_CellAttrMatName_Key);

    // Pre-compute axis dims here so the algorithm doesn't have to repeat it. Preflight has already
    // validated that binSizeDeg cleanly divides 360 and 180, so std::round is exact.
    inputValues.nphi1 = static_cast<int32_t>(std::round(360.0 / inputValues.binSizeDeg));
    inputValues.nPHI = static_cast<int32_t>(std::round(180.0 / inputValues.binSizeDeg));
    inputValues.nphi2 = inputValues.nphi1;
  }
  else
  {
    // Append mode: override target geometry / cell attr matrix name / bin size from the existing geometry.
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

  return ComputeODF(dataStructure, messageHandler, shouldCancel, &inputValues)();
}

//------------------------------------------------------------------------------
Result<Arguments> ComputeODFFilter::FromSIMPLJson(const nlohmann::json& json)
{
  Arguments args = ComputeODFFilter().getDefaultArguments();

  std::vector<Result<>> results;

  /* This is a NEW filter and has no SIMPL (DREAM3D v6) equivalent. */

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult), std::move(args));
}

} // namespace nx::core

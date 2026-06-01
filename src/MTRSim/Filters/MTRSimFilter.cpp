#include "MTRSimFilter.hpp"

#include "MTRSim/Filters/Algorithms/MTRSim.hpp"

#include "simplnx/Common/DataTypeUtilities.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/DynamicTableParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/MultiArraySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"
#include "simplnx/Parameters/util/DynamicTableInfo.hpp"

#include <fmt/format.h>

#include <chrono>
#include <cmath>
#include <random>
#include <vector>

using namespace nx::core;

namespace nx::core
{
//------------------------------------------------------------------------------
std::string MTRSimFilter::name() const
{
  return FilterTraits<MTRSimFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string MTRSimFilter::className() const
{
  return FilterTraits<MTRSimFilter>::className;
}

//------------------------------------------------------------------------------
Uuid MTRSimFilter::uuid() const
{
  return FilterTraits<MTRSimFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string MTRSimFilter::humanName() const
{
  return "Generate Synthetic Microtexture";
}

//------------------------------------------------------------------------------
std::vector<std::string> MTRSimFilter::defaultTags() const
{
  return {className(), "MTRSim", "Synthetic", "Microtexture", "Generate"};
}

//------------------------------------------------------------------------------
Parameters MTRSimFilter::parameters() const
{
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input ODF"});
  params.insert(std::make_unique<GeometrySelectionParameter>(k_InputOdfGeometry_Key, "Input ODF Geometry", "Image Geometry holding the ODF (from the Read/Compute ODF filters).", DataPath{},
                                                             GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<MultiArraySelectionParameter>(k_OdfComponentArrays_Key, "ODF Component Arrays", "Ordered list of per-component ODF cell arrays. Order maps to Volume Fraction columns.",
                                                              MultiArraySelectionParameter::ValueType{}, MultiArraySelectionParameter::AllowedTypes{IArray::ArrayType::DataArray}, GetAllNumericTypes(),
                                                              MultiArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insertSeparator(Parameters::Separator{"Simulation Parameters"});
  {
    DynamicTableInfo vfInfo;
    vfInfo.setRowsInfo(DynamicTableInfo::StaticVectorInfo(1));
    vfInfo.setColsInfo(DynamicTableInfo::DynamicVectorInfo(1, 3, "Comp {}"));
    params.insert(std::make_unique<DynamicTableParameter>(k_VolumeFractions_Key, "Volume Fraction", "One value per ODF component; must match the component count and sum to 1.0.", vfInfo));
  }
  {
    DynamicTableInfo thetaInfo;
    thetaInfo.setRowsInfo(DynamicTableInfo::DynamicVectorInfo(1, 2, "Gaussian {}"));
    thetaInfo.setColsInfo(DynamicTableInfo::StaticVectorInfo(DynamicTableInfo::HeadersListType{"theta_x", "theta_y", "theta_z"}));
    params.insert(std::make_unique<DynamicTableParameter>(k_ThetaList_Key, "Theta List",
                                                          "Correlation lengths [theta_x, theta_y, theta_z] per latent Gaussian. Needs >= (components - 1) rows. Same length unit as Physical "
                                                          "Size/Spacing.",
                                                          thetaInfo));
  }
  params.insert(std::make_unique<VectorFloat32Parameter>(k_PhysicalSize_Key, "Physical Size (microns)", "Domain extent X, Y, Z in microns. Set Z = 0 to generate a single-layer (2D) microstructure.", std::vector<float32>{38.1f, 12.7f, 0.0f},
                                                         std::vector<std::string>{"X", "Y", "Z"}));
  params.insert(std::make_unique<VectorFloat32Parameter>(k_PhysicalSpacing_Key, "Physical Spacing (microns)", "Voxel spacing X,Y,Z.", std::vector<float32>{0.02f, 0.02f, 0.02f},
                                                         std::vector<std::string>{"X", "Y", "Z"}));

  params.insertSeparator(Parameters::Separator{"Random Number Seed Parameters"});
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_UseSeed_Key, "Use Seed for Random Generation", "When true the user can supply a fixed seed.", false));
  params.insert(std::make_unique<NumberParameter<uint64>>(k_SeedValue_Key, "Seed Value", "The seed fed into the random generator.", std::mt19937::default_seed));
  params.insert(std::make_unique<DataObjectNameParameter>(k_SeedArrayName_Key, "Stored Seed Value Array Name", "Top-level array recording the seed used.", "MTRSim SeedValue"));

  params.insertSeparator(Parameters::Separator{"Output Data Object(s)"});
  params.insertLinkableParameter(
      std::make_unique<BoolParameter>(k_GeneratePolarColoring_Key, "Generate Polar Coloring", "Create a 3-component UInt8 RGB array using the MATLAB polar color mapping.", false));
  params.insert(std::make_unique<DataGroupCreationParameter>(k_OutputGeometry_Key, "Output Image Geometry", "Path of the new microstructure Image Geometry.", DataPath({"MTR Microstructure"})));
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellAttrMatName_Key, "Cell Attribute Matrix Name", "Name of the created cell AttributeMatrix.", "Cell Data"));
  params.insert(std::make_unique<DataObjectNameParameter>(k_MtrIdsArrayName_Key, "MTR Ids Array Name", "Int32 per-voxel MTR component id (1-based).", "MTRIds"));
  params.insert(std::make_unique<DataObjectNameParameter>(k_EulersArrayName_Key, "Euler Angles Array Name", "Float32 3-component Bunge Euler angles [rad].", "Eulers"));
  params.insert(std::make_unique<DataObjectNameParameter>(k_PolarColorsArrayName_Key, "Polar Colors Array Name", "UInt8 3-component RGB polar coloring.", "Polar Colors"));

  params.linkParameters(k_UseSeed_Key, k_SeedValue_Key, true);
  params.linkParameters(k_GeneratePolarColoring_Key, k_PolarColorsArrayName_Key, true);

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType MTRSimFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer MTRSimFilter::clone() const
{
  return std::make_unique<MTRSimFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult MTRSimFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                                     const ExecutionContext& executionContext) const
{
  auto pOdfArrays = filterArgs.value<MultiArraySelectionParameter::ValueType>(k_OdfComponentArrays_Key);
  auto pVolumeFractions = filterArgs.value<DynamicTableParameter::ValueType>(k_VolumeFractions_Key);
  auto pThetaList = filterArgs.value<DynamicTableParameter::ValueType>(k_ThetaList_Key);
  auto pSize = filterArgs.value<std::vector<float32>>(k_PhysicalSize_Key);
  auto pSpacing = filterArgs.value<std::vector<float32>>(k_PhysicalSpacing_Key);
  auto pGenPolar = filterArgs.value<bool>(k_GeneratePolarColoring_Key);
  auto pOutGeomPath = filterArgs.value<DataPath>(k_OutputGeometry_Key);
  auto pCellAttrMatName = filterArgs.value<std::string>(k_CellAttrMatName_Key);
  auto pMtrIdsName = filterArgs.value<std::string>(k_MtrIdsArrayName_Key);
  auto pEulersName = filterArgs.value<std::string>(k_EulersArrayName_Key);
  auto pPolarName = filterArgs.value<std::string>(k_PolarColorsArrayName_Key);
  auto pSeedArrayName = filterArgs.value<std::string>(k_SeedArrayName_Key);

  nx::core::Result<OutputActions> resultOutputActions;
  std::vector<PreflightValue> preflightUpdatedValues;

  const usize numComponents = pOdfArrays.size();
  if(numComponents < 2)
  {
    return {MakeErrorResult<OutputActions>(-13001, "MTRSim requires at least 2 ODF component arrays.")};
  }
  if(pVolumeFractions.size() != 1 || pVolumeFractions[0].size() != numComponents)
  {
    return {MakeErrorResult<OutputActions>(-13002, fmt::format("Volume Fraction must be 1 row x {} columns (one per ODF component).", numComponents))};
  }
  for(double v : pVolumeFractions[0])
  {
    if(v < 0.0 || v > 1.0)
    {
      return {MakeErrorResult<OutputActions>(-13007, fmt::format("Each Volume Fraction value must be in the range [0, 1] (got {:.4f}).", v))};
    }
  }
  double vfSum = 0.0;
  for(double v : pVolumeFractions[0])
  {
    vfSum += v;
  }
  if(std::abs(vfSum - 1.0) > 1.0e-3)
  {
    return {MakeErrorResult<OutputActions>(-13003, fmt::format("Volume Fraction values must sum to 1.0 (got {:.4f}).", vfSum))};
  }
  if(pThetaList.size() < numComponents - 1)
  {
    return {MakeErrorResult<OutputActions>(-13004, fmt::format("Theta List needs at least {} rows (components - 1).", numComponents - 1))};
  }
  for(const auto& row : pThetaList)
  {
    if(row.size() != 3)
    {
      return {MakeErrorResult<OutputActions>(-13005, "Each Theta List row must have exactly 3 columns.")};
    }
  }

  if(pSpacing[0] <= 0.0f || pSpacing[1] <= 0.0f)
  {
    return {MakeErrorResult<OutputActions>(-13006, "Physical Spacing X and Y must be greater than 0.")};
  }

  const auto dim = [](float len, float sp) { return static_cast<usize>(std::max(std::lround(len / sp), 1L)); };
  const usize nx = dim(pSize[0], pSpacing[0]);
  const usize ny = dim(pSize[1], pSpacing[1]);
  const usize nz = (pSize[2] <= 0.0f) ? 1 : dim(pSize[2], pSpacing[2]);

  const std::vector<usize> imageGeomDimsXYZ = {nx, ny, nz};
  const std::vector<float32> origin = {0.0f, 0.0f, 0.0f};
  const std::vector<float32> spacingXYZ = {pSpacing[0], pSpacing[1], pSpacing[2]};
  const std::vector<usize> tupleShapeZYX = {nz, ny, nx};

  resultOutputActions.value().appendAction(std::make_unique<CreateImageGeometryAction>(pOutGeomPath, imageGeomDimsXYZ, origin, spacingXYZ, pCellAttrMatName));

  const DataPath cellAttrMatPath = pOutGeomPath.createChildPath(pCellAttrMatName);
  resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::int32, tupleShapeZYX, std::vector<usize>{1}, cellAttrMatPath.createChildPath(pMtrIdsName)));
  resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::float32, tupleShapeZYX, std::vector<usize>{3}, cellAttrMatPath.createChildPath(pEulersName)));
  if(pGenPolar)
  {
    resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::uint8, tupleShapeZYX, std::vector<usize>{3}, cellAttrMatPath.createChildPath(pPolarName)));
  }
  resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::uint64, std::vector<usize>{1}, std::vector<usize>{1}, DataPath({pSeedArrayName})));

  preflightUpdatedValues.push_back({"Output Grid (X, Y, Z)", fmt::format("{} x {} x {}", nx, ny, nz)});
  preflightUpdatedValues.push_back({"Number of ODF Components", std::to_string(numComponents)});

  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> MTRSimFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                   const ExecutionContext& executionContext) const
{
  MTRSimInputValues inputValues;
  inputValues.inputOdfGeometryPath = filterArgs.value<DataPath>(k_InputOdfGeometry_Key);
  inputValues.odfComponentPaths = filterArgs.value<MultiArraySelectionParameter::ValueType>(k_OdfComponentArrays_Key);
  inputValues.volumeFractions = filterArgs.value<DynamicTableParameter::ValueType>(k_VolumeFractions_Key);
  inputValues.thetaList = filterArgs.value<DynamicTableParameter::ValueType>(k_ThetaList_Key);
  inputValues.physicalSize = filterArgs.value<std::vector<float32>>(k_PhysicalSize_Key);
  inputValues.physicalSpacing = filterArgs.value<std::vector<float32>>(k_PhysicalSpacing_Key);
  inputValues.generatePolarColoring = filterArgs.value<bool>(k_GeneratePolarColoring_Key);
  inputValues.outputGeometryPath = filterArgs.value<DataPath>(k_OutputGeometry_Key);
  inputValues.cellAttrMatName = filterArgs.value<std::string>(k_CellAttrMatName_Key);
  inputValues.mtrIdsArrayName = filterArgs.value<std::string>(k_MtrIdsArrayName_Key);
  inputValues.eulersArrayName = filterArgs.value<std::string>(k_EulersArrayName_Key);
  inputValues.polarColorsArrayName = filterArgs.value<std::string>(k_PolarColorsArrayName_Key);

  uint64 seed = filterArgs.value<uint64>(k_SeedValue_Key);
  if(!filterArgs.value<bool>(k_UseSeed_Key))
  {
    seed = static_cast<uint64>(std::chrono::steady_clock::now().time_since_epoch().count());
  }
  dataStructure.getDataRefAs<UInt64Array>(DataPath({filterArgs.value<std::string>(k_SeedArrayName_Key)}))[0] = seed;
  inputValues.seed = seed;

  return MTRSim(dataStructure, messageHandler, shouldCancel, &inputValues)();
}

//------------------------------------------------------------------------------
Result<Arguments> MTRSimFilter::FromSIMPLJson(const nlohmann::json& json)
{
  Arguments args = MTRSimFilter().getDefaultArguments();

  std::vector<Result<>> results;

  /* This is a NEW filter and has no SIMPL (DREAM3D v6) equivalent. */

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult), std::move(args));
}

} // namespace nx::core

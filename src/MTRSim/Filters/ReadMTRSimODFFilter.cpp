#include "ReadMTRSimODFFilter.hpp"

#include "MTRSim/Filters/Algorithms/ReadMTRSimODF.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"

#include "LibMTRSim/ODFFileIO.hpp"

#include <fmt/format.h>

#include <exception>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

using namespace nx::core;

namespace nx::core {
//------------------------------------------------------------------------------
std::string ReadMTRSimODFFilter::name() const {
  return FilterTraits<ReadMTRSimODFFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string ReadMTRSimODFFilter::className() const {
  return FilterTraits<ReadMTRSimODFFilter>::className;
}

//------------------------------------------------------------------------------
Uuid ReadMTRSimODFFilter::uuid() const {
  return FilterTraits<ReadMTRSimODFFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string ReadMTRSimODFFilter::humanName() const {
  return "Read MTRSim ODF (HDF5)";
}

//------------------------------------------------------------------------------
std::vector<std::string> ReadMTRSimODFFilter::defaultTags() const {
  return {className(), "IO", "Input", "Read", "MTRSim", "ODF"};
}

//------------------------------------------------------------------------------
Parameters ReadMTRSimODFFilter::parameters() const {
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<FileSystemPathParameter>(
      k_InputFile_Key, "Input ODF File (HDF5)",
      "MATLAB-format MTRSim ODF HDF5 file to read.", fs::path(""),
      FileSystemPathParameter::ExtensionsType{".h5", ".hdf5"},
      FileSystemPathParameter::PathType::InputFile));
  params.insert(std::make_unique<StringParameter>(
      k_Hdf5PathPrefix_Key, "HDF5 Path Prefix",
      "Internal HDF5 group path where the ODF data lives. MATLAB-generated "
      "files "
      "typically use '/ODF_best' (the default). If your file was saved in "
      "MATLAB "
      "under a variable named 'my_sample_ODF', set this to '/my_sample_ODF'.",
      "/ODF_best"));

  params.insertSeparator(Parameters::Separator{"Output Data Object(s)"});
  params.insert(std::make_unique<DataGroupCreationParameter>(
      k_OutputImageGeometry_Key, "Created Image Geometry",
      "Path at which the ImageGeom holding the ODF data will be created. One "
      "Float64 cell-data array is "
      "created per ODF component in the file.",
      DataPath({"ODF"})));
  params.insert(std::make_unique<DataObjectNameParameter>(
      k_CellAttrMatName_Key, "Cell Attribute Matrix Name",
      "Name of the cell AttributeMatrix created under the output ImageGeom; "
      "per-component ODF arrays are placed inside it.",
      "Cell Data"));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ReadMTRSimODFFilter::parametersVersion() const {
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ReadMTRSimODFFilter::clone() const {
  return std::make_unique<ReadMTRSimODFFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ReadMTRSimODFFilter::preflightImpl(
    const DataStructure &dataStructure, const Arguments &filterArgs,
    const MessageHandler &messageHandler, const std::atomic_bool &shouldCancel,
    const ExecutionContext &executionContext) const {
  auto pInputFile =
      filterArgs.value<FileSystemPathParameter::ValueType>(k_InputFile_Key);
  auto pHdf5PathPrefix =
      filterArgs.value<StringParameter::ValueType>(k_Hdf5PathPrefix_Key);
  auto pOutputImageGeomPath =
      filterArgs.value<DataPath>(k_OutputImageGeometry_Key);
  auto pCellAttrMatName = filterArgs.value<DataObjectNameParameter::ValueType>(
      k_CellAttrMatName_Key);

  nx::core::Result<OutputActions> resultOutputActions;
  std::vector<PreflightValue> preflightUpdatedValues;

  // Read and validate ODF metadata from the file. Any error from
  // readODFMetadata() is surfaced as a preflight failure so the UI shows it
  // before execute() is called.
  mtrsim::ODFFileMetadata metadata{};
  try {
    metadata = mtrsim::readODFMetadata(pInputFile, pHdf5PathPrefix);
  } catch (const std::exception &e) {
    return {MakeErrorResult<OutputActions>(
        -12001, fmt::format("ODF file read failed: {}", e.what()))};
  }

  // Axis mapping: phi1 is slowest-varying in the file, phi2 fastest.
  // ImageGeom XYZ must have phi2 as X (fastest), PHI as Y, phi1 as Z (slowest).
  const std::vector<usize> imageGeomDimsXYZ = {
      static_cast<usize>(metadata.dimsPhi1PHIPhi2[2]),
      static_cast<usize>(metadata.dimsPhi1PHIPhi2[1]),
      static_cast<usize>(metadata.dimsPhi1PHIPhi2[0])};
  const std::vector<float32> imageGeomOrigin = {0.0f, 0.0f, 0.0f};
  const std::vector<float32> imageGeomSpacingXYZ = {
      static_cast<float32>(metadata.spacingDegPhi1PHIPhi2[2]),
      static_cast<float32>(metadata.spacingDegPhi1PHIPhi2[1]),
      static_cast<float32>(metadata.spacingDegPhi1PHIPhi2[0])};

  // Tuple shape for the cell-data arrays is ZYX-ordered (simplnx convention).
  // This ordering also matches the row-major layout of
  // /ODF_best/component_*/ODFval so execute() can copy values straight across
  // without reshaping.
  const std::vector<usize> arrayTupleShapeZYX = {
      imageGeomDimsXYZ[2], imageGeomDimsXYZ[1], imageGeomDimsXYZ[0]};
  const std::vector<usize> componentShape = {1};

  {
    auto createImageGeomAction = std::make_unique<CreateImageGeometryAction>(
        pOutputImageGeomPath, imageGeomDimsXYZ, imageGeomOrigin,
        imageGeomSpacingXYZ, pCellAttrMatName);
    resultOutputActions.value().appendAction(std::move(createImageGeomAction));
  }

  const DataPath cellAttrMatPath =
      pOutputImageGeomPath.createChildPath(pCellAttrMatName);
  for (int64_t c = 0; c < metadata.numComponents; ++c) {
    DataPath componentPath =
        cellAttrMatPath.createChildPath(fmt::format("component_{}", c));
    auto createArrayAction = std::make_unique<CreateArrayAction>(
        DataType::float64, arrayTupleShapeZYX, componentShape, componentPath);
    resultOutputActions.value().appendAction(std::move(createArrayAction));
  }

  preflightUpdatedValues.push_back({"HDF5 Path Prefix", pHdf5PathPrefix});
  preflightUpdatedValues.push_back(
      {"Components Found", std::to_string(metadata.numComponents)});
  preflightUpdatedValues.push_back(
      {"Bins (phi1 x PHI x phi2)",
       fmt::format("{} x {} x {}", metadata.dimsPhi1PHIPhi2[0],
                   metadata.dimsPhi1PHIPhi2[1], metadata.dimsPhi1PHIPhi2[2])});
  preflightUpdatedValues.push_back(
      {"ImageGeom Dimensions (X, Y, Z)",
       fmt::format("{} x {} x {}", metadata.dimsPhi1PHIPhi2[2],
                   metadata.dimsPhi1PHIPhi2[1], metadata.dimsPhi1PHIPhi2[0])});
  preflightUpdatedValues.push_back(
      {"Spacing (X, Y, Z) [deg]",
       fmt::format("{:.4f}, {:.4f}, {:.4f}", metadata.spacingDegPhi1PHIPhi2[2],
                   metadata.spacingDegPhi1PHIPhi2[1],
                   metadata.spacingDegPhi1PHIPhi2[0])});

  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> ReadMTRSimODFFilter::executeImpl(
    DataStructure &dataStructure, const Arguments &filterArgs,
    const PipelineFilter *pipelineNode, const MessageHandler &messageHandler,
    const std::atomic_bool &shouldCancel,
    const ExecutionContext &executionContext) const {
  ReadMTRSimODFInputValues inputValues;

  inputValues.inputFile =
      filterArgs.value<FileSystemPathParameter::ValueType>(k_InputFile_Key);
  inputValues.hdf5PathPrefix =
      filterArgs.value<StringParameter::ValueType>(k_Hdf5PathPrefix_Key);
  inputValues.outputImageGeometryPath =
      filterArgs.value<DataPath>(k_OutputImageGeometry_Key);
  inputValues.cellAttrMatName =
      filterArgs.value<DataObjectNameParameter::ValueType>(
          k_CellAttrMatName_Key);

  return ReadMTRSimODF(dataStructure, messageHandler, shouldCancel,
                       &inputValues)();
}

//------------------------------------------------------------------------------
Result<Arguments>
ReadMTRSimODFFilter::FromSIMPLJson(const nlohmann::json &json) {
  Arguments args = ReadMTRSimODFFilter().getDefaultArguments();

  std::vector<Result<>> results;

  /* This is a NEW filter and has no SIMPL (DREAM3D v6) equivalent. */

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult),
                                    std::move(args));
}

} // namespace nx::core

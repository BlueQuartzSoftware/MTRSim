#include "WriteMTRSimODFFilter.hpp"

#include "MTRSim/Filters/Algorithms/WriteMTRSimODF.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/MultiArraySelectionParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"

#include <fmt/format.h>

#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

using namespace nx::core;

namespace nx::core {
//------------------------------------------------------------------------------
std::string WriteMTRSimODFFilter::name() const {
  return FilterTraits<WriteMTRSimODFFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string WriteMTRSimODFFilter::className() const {
  return FilterTraits<WriteMTRSimODFFilter>::className;
}

//------------------------------------------------------------------------------
Uuid WriteMTRSimODFFilter::uuid() const {
  return FilterTraits<WriteMTRSimODFFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string WriteMTRSimODFFilter::humanName() const {
  return "Write MTRSim ODF (HDF5)";
}

//------------------------------------------------------------------------------
std::vector<std::string> WriteMTRSimODFFilter::defaultTags() const {
  return {className(), "IO", "Output", "Write", "MTRSim", "ODF"};
}

//------------------------------------------------------------------------------
Parameters WriteMTRSimODFFilter::parameters() const {
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Output Parameter(s)"});
  params.insert(std::make_unique<FileSystemPathParameter>(
      k_OutputFile_Key, "Output ODF File (HDF5)",
      "Path of the MATLAB-format MTRSim ODF HDF5 file to create. An existing "
      "file is overwritten.",
      fs::path(""), FileSystemPathParameter::ExtensionsType{".h5", ".hdf5"},
      FileSystemPathParameter::PathType::OutputFile));
  params.insert(std::make_unique<StringParameter>(
      k_Hdf5PathPrefix_Key, "HDF5 Path Prefix",
      "Internal HDF5 group path where the ODF data will be written. Defaults "
      "to "
      "'/ODF_best' to match the MATLAB on-disk convention. Change this if you "
      "want "
      "the data saved under a different top-level group name (e.g. "
      "'/my_sample_ODF').",
      "/ODF_best"));

  params.insertSeparator(Parameters::Separator{"Input Data Objects"});
  params.insert(std::make_unique<GeometrySelectionParameter>(
      k_InputImageGeometry_Key, "Input Image Geometry",
      "ImageGeom whose (X,Y,Z) dimensions map to (phi2, PHI, phi1) on disk. "
      "The written axis ordering reverses this convention.",
      DataPath{},
      GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));

  params.insert(std::make_unique<MultiArraySelectionParameter>(
      k_ODFComponents_Key, "ODF Component Arrays",
      "Float64 single-component cell-data arrays on the input ImageGeom. Each "
      "array becomes one ODF component in the output file.",
      MultiArraySelectionParameter::ValueType{},
      MultiArraySelectionParameter::AllowedTypes{IArray::ArrayType::DataArray},
      MultiArraySelectionParameter::AllowedDataTypes{DataType::float64},
      MultiArraySelectionParameter::AllowedComponentShapes{{1}}));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType WriteMTRSimODFFilter::parametersVersion() const {
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer WriteMTRSimODFFilter::clone() const {
  return std::make_unique<WriteMTRSimODFFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult WriteMTRSimODFFilter::preflightImpl(
    const DataStructure &dataStructure, const Arguments &filterArgs,
    const MessageHandler &messageHandler, const std::atomic_bool &shouldCancel,
    const ExecutionContext &executionContext) const {
  auto pOutputFile =
      filterArgs.value<FileSystemPathParameter::ValueType>(k_OutputFile_Key);
  auto pHdf5PathPrefix =
      filterArgs.value<StringParameter::ValueType>(k_Hdf5PathPrefix_Key);
  auto pInputImageGeomPath =
      filterArgs.value<DataPath>(k_InputImageGeometry_Key);
  auto pODFComponents =
      filterArgs.value<MultiArraySelectionParameter::ValueType>(
          k_ODFComponents_Key);

  nx::core::Result<OutputActions> resultOutputActions;
  std::vector<PreflightValue> preflightUpdatedValues;

  if (pODFComponents.empty()) {
    return {MakeErrorResult<OutputActions>(
        -12100, "Select at least one ODF component array.")};
  }

  const auto *geom = dataStructure.getDataAs<ImageGeom>(pInputImageGeomPath);
  if (geom == nullptr) {
    return {MakeErrorResult<OutputActions>(
        -12101, "Input geometry must be an ImageGeom.")};
  }

  const usize numX = geom->getNumXCells();
  const usize numY = geom->getNumYCells();
  const usize numZ = geom->getNumZCells();
  if (numX < 2 || numY < 2 || numZ < 2) {
    return {MakeErrorResult<OutputActions>(
        -12105, "Degenerate geometry rejected (< 2 bins on some axis).")};
  }

  const usize expectedTuples = numX * numY * numZ;
  const auto geomPathVec = pInputImageGeomPath.getPathVector();

  for (const auto &arrayPath : pODFComponents) {
    // Verify the array lives under the selected ImageGeom (its path must start
    // with geomPath).
    const auto arrayPathVec = arrayPath.getPathVector();
    bool onGeom = (arrayPathVec.size() > geomPathVec.size());
    if (onGeom) {
      for (usize i = 0; i < geomPathVec.size(); ++i) {
        if (arrayPathVec[i] != geomPathVec[i]) {
          onGeom = false;
          break;
        }
      }
    }
    if (!onGeom) {
      return {MakeErrorResult<OutputActions>(
          -12102,
          fmt::format("Selected array '{}' does not belong to the input "
                      "ImageGeom '{}'.",
                      arrayPath.toString(), pInputImageGeomPath.toString()))};
    }

    const auto *arr = dataStructure.getDataAs<Float64Array>(arrayPath);
    if (arr == nullptr) {
      return {MakeErrorResult<OutputActions>(
          -12103,
          fmt::format("Selected array '{}' must be a Float64 DataArray.",
                      arrayPath.toString()))};
    }

    if (arr->getNumberOfTuples() != expectedTuples) {
      return {MakeErrorResult<OutputActions>(
          -12104, fmt::format("Selected array '{}' has {} tuples but the "
                              "ImageGeom has {} cells.",
                              arrayPath.toString(), arr->getNumberOfTuples(),
                              expectedTuples))};
    }
  }

  // Populate a preview of what execute() will write. We can't probe a file that
  // doesn't exist yet, so every field is derived from the selected geometry and
  // the array list.
  const auto spacing = geom->getSpacing();
  preflightUpdatedValues.push_back({"HDF5 Path Prefix", pHdf5PathPrefix});
  preflightUpdatedValues.push_back(
      {"Components to Write", std::to_string(pODFComponents.size())});
  preflightUpdatedValues.push_back(
      {"Tuples per Component", std::to_string(expectedTuples)});
  preflightUpdatedValues.push_back(
      {"Grid (phi1 x PHI x phi2)",
       fmt::format("{} x {} x {}", numZ, numY, numX)});
  preflightUpdatedValues.push_back(
      {"Spacing (phi1, PHI, phi2) [deg]",
       fmt::format("{:.4f}, {:.4f}, {:.4f}", spacing[2], spacing[1],
                   spacing[0])});

  // Estimate output size: N components * tuples * 8 bytes + some overhead for
  // bin arrays
  const std::size_t estBytes =
      static_cast<std::size_t>(pODFComponents.size()) *
          static_cast<std::size_t>(expectedTuples) * sizeof(double) +
      static_cast<std::size_t>(numX + numY + numZ + 3) * sizeof(double);
  preflightUpdatedValues.push_back(
      {"Estimated File Size [bytes]", std::to_string(estBytes)});

  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> WriteMTRSimODFFilter::executeImpl(
    DataStructure &dataStructure, const Arguments &filterArgs,
    const PipelineFilter *pipelineNode, const MessageHandler &messageHandler,
    const std::atomic_bool &shouldCancel,
    const ExecutionContext &executionContext) const {
  WriteMTRSimODFInputValues inputValues;

  inputValues.outputFile =
      filterArgs.value<FileSystemPathParameter::ValueType>(k_OutputFile_Key);
  inputValues.hdf5PathPrefix =
      filterArgs.value<StringParameter::ValueType>(k_Hdf5PathPrefix_Key);
  inputValues.inputImageGeometry =
      filterArgs.value<DataPath>(k_InputImageGeometry_Key);
  inputValues.odfComponents =
      filterArgs.value<MultiArraySelectionParameter::ValueType>(
          k_ODFComponents_Key);

  return WriteMTRSimODF(dataStructure, messageHandler, shouldCancel,
                        &inputValues)();
}

//------------------------------------------------------------------------------
Result<Arguments>
WriteMTRSimODFFilter::FromSIMPLJson(const nlohmann::json &json) {
  Arguments args = WriteMTRSimODFFilter().getDefaultArguments();

  std::vector<Result<>> results;

  /* This is a NEW filter and has no SIMPL (DREAM3D v6) equivalent. */

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult),
                                    std::move(args));
}

} // namespace nx::core

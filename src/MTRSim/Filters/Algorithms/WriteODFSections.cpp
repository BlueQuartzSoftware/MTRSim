#include "WriteODFSections.hpp"

#include "EbsdLib/Utilities/ODFSectionCompositor.h"
#include "EbsdLib/Utilities/PngWriter.h"

#include "simplnx/Common/AtomicFile.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/StringArray.hpp"
#include "simplnx/Utilities/ColorTableUtilities.hpp"
#include "simplnx/Utilities/FilterUtilities.hpp"

#include <fmt/format.h>

#include <algorithm>

using namespace nx::core;

namespace
{
struct CopyODFToDouble
{
  template <class T>
  void operator()(const IDataArray& input, ebsdlib::DoubleArrayType& output, const std::atomic_bool& shouldCancel) const
  {
    const auto& array = dynamic_cast<const DataArray<T>&>(input);
    // Bounded batches allow cancellation without an atomic load for each value.
    constexpr usize k_BatchSize = 65536;
    for(usize begin = 0; begin < array.getSize(); begin += k_BatchSize)
    {
      if(shouldCancel)
      {
        return;
      }
      const usize end = begin + std::min(k_BatchSize, array.getSize() - begin);
      for(usize index = begin; index < end; ++index)
      {
        output[index] = static_cast<double>(array[index]);
      }
    }
  }
};
} // namespace

WriteODFSections::WriteODFSections(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel, const WriteODFSectionsInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_MessageHandler(messageHandler)
, m_ShouldCancel(shouldCancel)
, m_InputValues(inputValues)
{
}

WriteODFSections::~WriteODFSections() noexcept = default;

Result<> WriteODFSections::operator()()
{
  if(m_ShouldCancel)
  {
    return {};
  }
  // Synchronous phase messages let the handler request cancellation before the next stage starts.
  m_MessageHandler("Converting ODF values");
  if(m_ShouldCancel)
  {
    return {};
  }

  const auto& values = *m_InputValues;
  const auto& inputArray = m_DataStructure.getDataRefAs<IDataArray>(values.odfArrayPath);
  auto ebsdValues = ebsdlib::DoubleArrayType::CreateArray(inputArray.getNumberOfTuples(), {1ULL}, values.odfArrayPath.toString(), true);
  ExecuteDataFunction(CopyODFToDouble{}, inputArray.getDataType(), inputArray, *ebsdValues, m_ShouldCancel);
  if(m_ShouldCancel)
  {
    return {};
  }

  auto controls = ColorTableUtilities::ExtractControlPoints(values.colorTablePreset);
  if(controls.invalid())
  {
    return MakeErrorResult(-12322, fmt::format("Could not resolve color preset '{}' for ODF '{}' and phase {}: {}", values.colorTablePreset, values.odfArrayPath.toString(), values.phaseIndex,
                                               controls.errors().front().message));
  }
  const auto& imageGeom = m_DataStructure.getDataRefAs<ImageGeom>(values.inputImageGeometryPath);
  const auto& crystalStructures = m_DataStructure.getDataRefAs<UInt32Array>(values.crystalStructuresPath);
  const auto& materialNames = m_DataStructure.getDataRefAs<StringArray>(values.materialNamesPath);
  const auto origin = imageGeom.getOrigin();
  const auto spacing = imageGeom.getSpacing();
  ebsdlib::ODFSectionConfiguration config;
  config.grid.values = ebsdValues.get();
  // Image geometry XYZ corresponds to phi2, PHI, phi1; the borrowed values retain phi2-fastest storage.
  config.grid.dimensions = {imageGeom.getNumZCells(), imageGeom.getNumYCells(), imageGeom.getNumXCells()};
  config.grid.originDeg = {origin[2], origin[1], origin[0]};
  config.grid.spacingDeg = {spacing[2], spacing[1], spacing[0]};
  config.grid.units = values.inputUnits == 0 ? ebsdlib::ODFValueUnits::MUD : ebsdlib::ODFValueUnits::CountDensity;
  config.grid.laueOpsIndex = crystalStructures[static_cast<usize>(values.phaseIndex)];
  config.sectionCount = static_cast<usize>(values.numberOfSections);
  config.sectionsPerRow = static_cast<usize>(values.sectionsPerRow);
  config.sectionWidth = values.sectionWidth;
  config.scaleMode = values.scaleMode == 0 ? ebsdlib::ODFScaleMode::Automatic : ebsdlib::ODFScaleMode::Manual;
  config.manualMinimumMUD = values.manualMinimumMUD;
  config.manualMaximumMUD = values.manualMaximumMUD;
  // SIMPLNX presets can use arbitrary scalar positions; EbsdLib accepts normalized positions.
  const auto binPoints = ColorTableUtilities::NormalizeBinPoints(controls.value());
  config.colorControlPoints = std::move(controls.value());
  for(usize colorIndex = 0; colorIndex < binPoints.size(); ++colorIndex)
  {
    config.colorControlPoints[colorIndex * 4] = binPoints[colorIndex];
  }
  config.title = values.title;
  config.materialName = materialNames[static_cast<usize>(values.phaseIndex)];
  config.phaseNumber = values.phaseIndex;

  m_MessageHandler("Composing ODF sections");
  if(m_ShouldCancel)
  {
    return {};
  }
  const auto result = ebsdlib::ODFSectionCompositor{}.generateCompositeImage(config);
  if(!result)
  {
    return MakeErrorResult(-12321,
                           fmt::format("Could not render ODF '{}' for phase {}: {} (EbsdLib code {}).", values.odfArrayPath.toString(), values.phaseIndex, result.errorMessage, result.errorCode));
  }
  m_MessageHandler(fmt::format("ODF section scale: minimum MUD = {}, maximum MUD = {}.", result.appliedMinimumMUD, result.appliedMaximumMUD));
  if(m_ShouldCancel)
  {
    return {};
  }

  if(values.saveAsImageGeometry)
  {
    m_MessageHandler("Copying ODF section RGB pixels");
    if(m_ShouldCancel)
    {
      return {};
    }
    const auto rgbPath = values.outputImageGeometryPath.createChildPath("Cell Data").createChildPath("ODF Sections RGB");
    auto& outputRgb = m_DataStructure.getDataRefAs<UInt8Array>(rgbPath);
    const auto* rgba = result.image->getPointer(0);
    const usize width = static_cast<usize>(result.width);
    for(usize row = 0; row < static_cast<usize>(result.height); ++row)
    {
      if(m_ShouldCancel)
      {
        return {};
      }
      for(usize column = 0; column < width; ++column)
      {
        const usize tuple = row * width + column;
        outputRgb[tuple * 3] = rgba[tuple * 4];
        outputRgb[tuple * 3 + 1] = rgba[tuple * 4 + 1];
        outputRgb[tuple * 3 + 2] = rgba[tuple * 4 + 2];
      }
    }
  }
  if(values.writePng)
  {
    m_MessageHandler("Writing ODF section PNG");
    if(m_ShouldCancel)
    {
      return {};
    }
    auto atomicResult = AtomicFile::Create(values.outputPngPath);
    if(atomicResult.invalid())
    {
      return ConvertResult(std::move(atomicResult));
    }
    AtomicFile atomicFile = std::move(atomicResult.value());
    const auto pngResult = PngWriter::WriteColorImage(atomicFile.tempFilePath().string(), result.width, result.height, 4, result.image->getPointer(0));
    if(pngResult.first < 0)
    {
      return MakeErrorResult(-12330, fmt::format("Could not write ODF section PNG '{}': {} (writer code {}).", values.outputPngPath.string(), pngResult.second, pngResult.first));
    }
    if(m_ShouldCancel)
    {
      return {};
    }
    return atomicFile.commit();
  }
  return {};
}

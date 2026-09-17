#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include <filesystem>
#include <string>

namespace nx::core
{
/**
 * @struct WriteODFSectionsInputValues
 * @brief Contains the resolved values for ODF section rendering.
 */
struct MTRSIM_EXPORT WriteODFSectionsInputValues
{
  DataPath inputImageGeometryPath;
  DataPath odfArrayPath;
  uint64 inputUnits = 0;
  DataPath crystalStructuresPath;
  DataPath materialNamesPath;
  int32 phaseIndex = 1;
  std::string title;
  int32 numberOfSections = 6;
  int32 sectionsPerRow = 3;
  int32 sectionWidth = 512;
  std::string colorTablePreset = "Rainbow";
  uint64 scaleMode = 0;
  float64 manualMinimumMUD = 0.0;
  float64 manualMaximumMUD = 1.0;
  bool writePng = true;
  std::filesystem::path outputPngPath;
  bool saveAsImageGeometry = true;
  DataPath outputImageGeometryPath = DataPath({"ODF Sections"});
};

/**
 * @class WriteODFSections
 * @brief Renders a standard ODF section montage from the resolved filter inputs.
 */
class MTRSIM_EXPORT WriteODFSections
{
public:
  /**
   * @brief Constructs the ODF section rendering algorithm.
   * @param dataStructure Contains the source ODF and destination image.
   * @param messageHandler Receives execution messages.
   * @param shouldCancel Indicates that the algorithm must stop.
   * @param inputValues Contains the resolved filter inputs for the call lifetime.
   */
  WriteODFSections(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel, const WriteODFSectionsInputValues* inputValues);
  /**
   * @brief Releases the algorithm resources.
   */
  ~WriteODFSections() noexcept;

  WriteODFSections(const WriteODFSections&) = delete;
  WriteODFSections(WriteODFSections&&) noexcept = delete;
  WriteODFSections& operator=(const WriteODFSections&) = delete;
  WriteODFSections& operator=(WriteODFSections&&) noexcept = delete;

  /**
   * @brief Renders the configured outputs.
   * @return Success or a rendering error.
   */
  Result<> operator()();

private:
  DataStructure& m_DataStructure;
  const IFilter::MessageHandler& m_MessageHandler;
  const std::atomic_bool& m_ShouldCancel;
  const WriteODFSectionsInputValues* m_InputValues = nullptr;
};
} // namespace nx::core

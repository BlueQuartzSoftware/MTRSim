#pragma once

#include "EbsdLib/Utilities/ODFSectionCompositor.h"
#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include <filesystem>
#include <string>

namespace nx::core
{
/**
 * @struct GenerateODFSectionsInputValues
 * @brief Contains the resolved values for ODF section rendering.
 */
struct MTRSIM_EXPORT GenerateODFSectionsInputValues
{
  DataPath inputImageGeometryPath;
  DataPath odfArrayPath;
  ebsdlib::ODFValueUnits units = ebsdlib::ODFValueUnits::MUD;
  uint32 crystalStructure = 0;
  std::string materialName;
  int32 phaseIndex = 1;
  std::string title;
  int32 numberOfSections = 6;
  int32 sectionsPerRow = 3;
  int32 sectionWidth = 512;
  std::string colorTablePreset = "Rainbow Blended White";
  ebsdlib::ODFScaleMode scaleMode = ebsdlib::ODFScaleMode::Automatic;
  float64 manualMinimumMUD = 0.0;
  float64 manualMaximumMUD = 1.0;
  bool writePng = false;
  std::filesystem::path outputPngPath;
  bool saveAsImageGeometry = true;
  DataPath outputImageGeometryPath = DataPath({"ODF Sections"});
};

/**
 * @class GenerateODFSections
 * @brief Renders a standard ODF section montage from the resolved presentation inputs.
 */
class MTRSIM_EXPORT GenerateODFSections
{
public:
  /**
   * @brief Constructs the ODF section rendering algorithm.
   * @param dataStructure Contains the source ODF and destination image.
   * @param messageHandler Receives execution messages.
   * @param shouldCancel Indicates that the algorithm must stop.
   * @param inputValues Contains the resolved presentation inputs for the call lifetime.
   */
  GenerateODFSections(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel, const GenerateODFSectionsInputValues* inputValues);
  /**
   * @brief Releases the algorithm resources.
   */
  ~GenerateODFSections() noexcept;

  GenerateODFSections(const GenerateODFSections&) = delete;
  GenerateODFSections(GenerateODFSections&&) noexcept = delete;
  GenerateODFSections& operator=(const GenerateODFSections&) = delete;
  GenerateODFSections& operator=(GenerateODFSections&&) noexcept = delete;

  /**
   * @brief Renders the configured outputs.
   * @pre Preflight actions have created the optional RGB geometry.
   * @return Success, cancellation, or a contextual rendering or file error.
   * @note Cancellation can leave a partial RGB array. The PNG destination changes only after successful encoding and atomic commit.
   */
  Result<> operator()();

private:
  DataStructure& m_DataStructure;
  const IFilter::MessageHandler& m_MessageHandler;
  const std::atomic_bool& m_ShouldCancel;
  const GenerateODFSectionsInputValues* m_InputValues = nullptr;
};
} // namespace nx::core

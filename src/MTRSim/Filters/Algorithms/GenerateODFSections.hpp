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
  DataPath input_image_geometry_path;
  DataPath odf_array_path;
  ebsdlib::ODFValueUnits units = ebsdlib::ODFValueUnits::MUD;
  uint32 crystal_structure = 0;
  std::string material_name;
  int32 phase_index = 1;
  std::string title;
  int32 number_of_sections = 6;
  int32 sections_per_row = 3;
  int32 section_width = 512;
  std::string color_table_preset = "Rainbow Blended White";
  ebsdlib::ODFScaleMode scale_mode = ebsdlib::ODFScaleMode::Automatic;
  float64 manual_minimum_mud = 0.0;
  float64 manual_maximum_mud = 1.0;
  bool write_png = false;
  std::filesystem::path output_png_path;
  bool save_as_image_geometry = true;
  DataPath output_image_geometry_path = DataPath({"ODF Sections"});
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
  ~GenerateODFSections() noexcept = default;

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

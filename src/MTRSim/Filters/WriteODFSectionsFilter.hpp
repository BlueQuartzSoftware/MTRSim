#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/Common/StringLiteral.hpp"
#include "simplnx/Filter/FilterTraits.hpp"
#include "simplnx/Filter/IFilter.hpp"

namespace nx::core
{
/**
 * @class WriteODFSectionsFilter
 * @brief Creates a standard constant-phi2 ODF section montage from a full Bunge Euler grid.
 */
class MTRSIM_EXPORT WriteODFSectionsFilter : public IFilter
{
public:
  /**
   * @brief Constructs the filter.
   */
  WriteODFSectionsFilter() = default;
  /**
   * @brief Releases the filter resources.
   */
  ~WriteODFSectionsFilter() noexcept override = default;

  WriteODFSectionsFilter(const WriteODFSectionsFilter&) = delete;
  WriteODFSectionsFilter(WriteODFSectionsFilter&&) noexcept = delete;
  WriteODFSectionsFilter& operator=(const WriteODFSectionsFilter&) = delete;
  WriteODFSectionsFilter& operator=(WriteODFSectionsFilter&&) noexcept = delete;

  static constexpr StringLiteral k_InputImageGeometryPath_Key = "input_image_geometry_path";
  static constexpr StringLiteral k_ODFArrayPath_Key = "odf_array_path";
  static constexpr StringLiteral k_InputUnits_Key = "input_units_index";
  static constexpr StringLiteral k_CrystalStructuresPath_Key = "crystal_structures_path";
  static constexpr StringLiteral k_MaterialNamesPath_Key = "material_names_path";
  static constexpr StringLiteral k_PhaseIndex_Key = "phase_index";
  static constexpr StringLiteral k_Title_Key = "title";
  static constexpr StringLiteral k_NumSections_Key = "number_of_sections";
  static constexpr StringLiteral k_SectionsPerRow_Key = "sections_per_row";
  static constexpr StringLiteral k_SectionWidth_Key = "section_width";
  static constexpr StringLiteral k_ColorTablePreset_Key = "color_table_preset";
  static constexpr StringLiteral k_ScaleMode_Key = "scale_mode_index";
  static constexpr StringLiteral k_ManualMinimum_Key = "manual_minimum_mud";
  static constexpr StringLiteral k_ManualMaximum_Key = "manual_maximum_mud";
  static constexpr StringLiteral k_WritePng_Key = "write_png";
  static constexpr StringLiteral k_OutputPngPath_Key = "output_png_file";
  static constexpr StringLiteral k_SaveAsImageGeometry_Key = "save_as_image_geometry";
  static constexpr StringLiteral k_OutputImageGeometryPath_Key = "output_image_geometry_path";

  /**
   * @brief Returns the filter registry name.
   * @return Filter registry name.
   */
  std::string name() const override;
  /**
   * @brief Returns the filter class name.
   * @return Filter class name.
   */
  std::string className() const override;
  /**
   * @brief Returns the stable filter identifier.
   * @return Filter UUID.
   */
  Uuid uuid() const override;
  /**
   * @brief Returns the user-facing filter name.
   * @return User-facing filter name.
   */
  std::string humanName() const override;
  /**
   * @brief Returns the default discovery tags.
   * @return Tags for filter search and categorization.
   */
  std::vector<std::string> defaultTags() const override;
  /**
   * @brief Defines the ODF, presentation, and output parameters.
   * @return Filter parameter collection.
   */
  Parameters parameters() const override;
  /**
   * @brief Returns the parameter schema version.
   * @return Parameter schema version.
   */
  VersionType parametersVersion() const override;
  /**
   * @brief Creates a new filter instance.
   * @return New filter instance.
   */
  UniquePointer clone() const override;

protected:
  /**
   * @brief Validates the ODF grid and declares the selected montage outputs.
   * @param dataStructure Contains the input ODF and phase metadata.
   * @param filterArgs Contains the resolved filter arguments.
   * @param messageHandler Receives preflight messages.
   * @param shouldCancel Indicates that the filter must stop.
   * @param executionContext Resolves relative filesystem paths.
   * @return Output actions and ODF section preview values, or validation errors.
   */
  PreflightResult preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                const ExecutionContext& executionContext) const override;

  /**
   * @brief Dispatches ODF section rendering.
   * @param dataStructure Contains the input and created output objects.
   * @param filterArgs Contains the resolved filter arguments.
   * @param pipelineNode Identifies the active pipeline node.
   * @param messageHandler Receives execution messages.
   * @param shouldCancel Indicates that the filter must stop.
   * @param executionContext Resolves relative filesystem paths.
   * @return Success or an execution error.
   */
  Result<> executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                       const ExecutionContext& executionContext) const override;
};
} // namespace nx::core

SIMPLNX_DEF_FILTER_TRAITS(nx::core, WriteODFSectionsFilter, "5ec65380-b1a8-4aff-8bcd-310b5067b236");

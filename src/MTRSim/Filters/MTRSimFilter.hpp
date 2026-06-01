#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/Common/StringLiteral.hpp"
#include "simplnx/Filter/FilterTraits.hpp"
#include "simplnx/Filter/IFilter.hpp"

namespace nx::core {
/**
 * @class MTRSimFilter
 * @brief Generates a synthetic microtexture (MTR) microstructure from an input
 * ODF and a set of simulation parameters, producing a new ImageGeom with
 * per-voxel MTR Ids, Euler angles, and optional polar coloring.
 */
class MTRSIM_EXPORT MTRSimFilter : public IFilter {
public:
  MTRSimFilter() = default;
  ~MTRSimFilter() noexcept override = default;

  MTRSimFilter(const MTRSimFilter &) = delete;
  MTRSimFilter(MTRSimFilter &&) noexcept = delete;

  MTRSimFilter &operator=(const MTRSimFilter &) = delete;
  MTRSimFilter &operator=(MTRSimFilter &&) noexcept = delete;

  // Parameter Keys
  static inline constexpr StringLiteral k_InputOdfGeometry_Key =
      "input_odf_geometry";
  static inline constexpr StringLiteral k_OdfComponentArrays_Key =
      "odf_component_arrays";
  static inline constexpr StringLiteral k_VolumeFractions_Key =
      "volume_fractions";
  static inline constexpr StringLiteral k_ThetaList_Key = "theta_list";
  static inline constexpr StringLiteral k_PhysicalSize_Key = "physical_size";
  static inline constexpr StringLiteral k_PhysicalSpacing_Key =
      "physical_spacing";
  static inline constexpr StringLiteral k_UseSeed_Key = "use_seed";
  static inline constexpr StringLiteral k_SeedValue_Key = "seed_value";
  static inline constexpr StringLiteral k_SeedArrayName_Key = "seed_array_name";
  static inline constexpr StringLiteral k_GeneratePolarColoring_Key =
      "generate_polar_coloring";
  static inline constexpr StringLiteral k_OutputGeometry_Key =
      "output_geometry";
  static inline constexpr StringLiteral k_CellAttrMatName_Key =
      "cell_attribute_matrix_name";
  static inline constexpr StringLiteral k_MtrIdsArrayName_Key =
      "mtr_ids_array_name";
  static inline constexpr StringLiteral k_EulersArrayName_Key =
      "eulers_array_name";
  static inline constexpr StringLiteral k_PolarColorsArrayName_Key =
      "polar_colors_array_name";

  /**
   * @brief Reads SIMPL json and converts it simplnx Arguments.
   * @param json
   * @return Result<Arguments>
   */
  static Result<Arguments> FromSIMPLJson(const nlohmann::json &json);

  /**
   * @brief Returns the name of the filter.
   * @return
   */
  std::string name() const override;

  /**
   * @brief Returns the C++ classname of this filter.
   * @return
   */
  std::string className() const override;

  /**
   * @brief Returns the uuid of the filter.
   * @return
   */
  Uuid uuid() const override;

  /**
   * @brief Returns the human readable name of the filter.
   * @return
   */
  std::string humanName() const override;

  /**
   * @brief Returns the default tags for this filter.
   * @return
   */
  std::vector<std::string> defaultTags() const override;

  /**
   * @brief Returns the parameters of the filter (i.e. its inputs)
   * @return
   */
  Parameters parameters() const override;

  /**
   * @brief Returns parameters version integer.
   * Initial version should always be 1.
   * Should be incremented everytime the parameters change.
   * @return VersionType
   */
  VersionType parametersVersion() const override;

  /**
   * @brief Returns a copy of the filter.
   * @return
   */
  UniquePointer clone() const override;

protected:
  /**
   * @brief Takes in a DataStructure and checks that the filter can be run on it
   * with the given arguments. Returns any warnings/errors. Also returns the
   * changes that would be applied to the DataStructure. Some parts of the
   * actions may not be completely filled out if all the required information is
   * not available at preflight time.
   * @param dataStructure The input DataStructure instance
   * @param filterArgs These are the input values for each parameter that is
   * required for the filter
   * @param messageHandler The MessageHandler object
   * @param shouldCancel Atomic boolean value that can be checked to cancel the
   * filter
   * @param executionContext The ExecutionContext that can be used to determine
   * the correct absolute path from a relative path
   * @return Returns a Result object with error or warning values if any of
   * those occurred during execution of this function
   */
  PreflightResult
  preflightImpl(const DataStructure &dataStructure, const Arguments &filterArgs,
                const MessageHandler &messageHandler,
                const std::atomic_bool &shouldCancel,
                const ExecutionContext &executionContext) const override;

  /**
   * @brief Applies the filter's algorithm to the DataStructure with the given
   * arguments. Returns any warnings/errors. On failure, there is no guarantee
   * that the DataStructure is in a correct state.
   * @param dataStructure The input DataStructure instance
   * @param filterArgs These are the input values for each parameter that is
   * required for the filter
   * @param pipelineNode The node in the pipeline that is being executed
   * @param messageHandler The MessageHandler object
   * @param shouldCancel Atomic boolean value that can be checked to cancel the
   * filter
   * @param executionContext The ExecutionContext that can be used to determine
   * the correct absolute path from a relative path
   * @return Returns a Result object with error or warning values if any of
   * those occurred during execution of this function
   */
  Result<> executeImpl(DataStructure &dataStructure,
                       const Arguments &filterArgs,
                       const PipelineFilter *pipelineNode,
                       const MessageHandler &messageHandler,
                       const std::atomic_bool &shouldCancel,
                       const ExecutionContext &executionContext) const override;
};
} // namespace nx::core

SIMPLNX_DEF_FILTER_TRAITS(nx::core, MTRSimFilter,
                          "f7f7a330-4bff-4a42-a573-09117a89a0a0");

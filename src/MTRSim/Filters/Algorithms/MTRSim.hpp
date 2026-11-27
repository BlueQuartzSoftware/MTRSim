#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include "LibMTRSim/MTRSimDriver.hpp"

#include <vector>

namespace nx::core
{

struct MTRSIM_EXPORT MTRSimInputValues
{
  DataPath input_odf_geometry_path;
  std::vector<DataPath> odf_component_paths;
  std::vector<std::vector<double>> volume_fractions; // 1 row x N cols
  std::vector<std::vector<double>> theta_list;       // M rows x 3 cols
  std::vector<float32> physical_size;                // [x,y,z] microns
  std::vector<float32> physical_spacing;             // [x,y,z] microns
  uint64 seed{};
  bool generate_polar_coloring{};
  DataPath output_geometry_path;
  std::string cell_attr_mat_name;
  std::string mtr_ids_array_name;
  std::string eulers_array_name;
  std::string polar_colors_array_name;
};

/**
 * @class MTRSim
 * @brief Algorithm that generates a synthetic microtexture (MTR) microstructure
 * from an input ODF and a set of simulation parameters. The output ImageGeom
 * and its cell arrays are created by the filter's preflight; this algorithm
 * fills them.
 */
class MTRSIM_EXPORT MTRSim
{
public:
  MTRSim(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, MTRSimInputValues* inputValues);
  ~MTRSim() noexcept = default;

  MTRSim(const MTRSim&) = delete;
  MTRSim(MTRSim&&) noexcept = delete;
  MTRSim& operator=(const MTRSim&) = delete;
  MTRSim& operator=(MTRSim&&) noexcept = delete;

  Result<> operator()();

private:
  Result<> applyPolarColoring(const mtrsim::MTRSimResult& sim, const DataPath& cellAttrMatPath);
  DataStructure& m_DataStructure;
  const MTRSimInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core

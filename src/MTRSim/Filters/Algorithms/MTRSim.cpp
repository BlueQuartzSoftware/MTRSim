#include "MTRSim.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"

#include "LibMTRSim/IPFMapper.hpp"
#include "LibMTRSim/ISimulationObserver.hpp"
#include "LibMTRSim/MTRSimDriver.hpp"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <atomic>
#include <cstdint>
#include <exception>
#include <random>
#include <vector>

using namespace nx::core;

namespace
{
/**
 * @brief Adapts mtrsim::ISimulationObserver to the simplnx filter message
 * handler and cancel flag, so simulateMTR can report progress and be cancelled.
 */
class FilterObserver : public mtrsim::ISimulationObserver
{
public:
  FilterObserver(const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel)
  : m_MessageHandler(messageHandler)
  , m_ShouldCancel(shouldCancel)
  {
  }

  void updateProgress(int64_t done, int64_t total, const std::string& message) override
  {
    const int32 progress = (total > 0) ? static_cast<int32>(done * 100 / total) : 0;
    m_MessageHandler.sendProgressMessage(message, progress);
  }

  void info(const std::string& message) override
  {
    m_MessageHandler.sendInfoMessage(message);
  }

  [[nodiscard]] bool shouldCancel() const override
  {
    return m_ShouldCancel.load();
  }

private:
  const IFilter::MessageHandler& m_MessageHandler;
  const std::atomic_bool& m_ShouldCancel;
};
} // namespace

// -----------------------------------------------------------------------------
MTRSim::MTRSim(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, MTRSimInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Result<> MTRSim::operator()()
{
  m_MessageHandler.sendInfoMessage("Reading ODF geometry...");

  // 1. ODF grid geometry (degrees) from the input ODF ImageGeom.
  const auto& odfGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->input_odf_geometry_path);
  const SizeVec3Type odfDims = odfGeom.getDimensions();  // X=phi2, Y=PHI, Z=phi1
  const FloatVec3Type odfSpacing = odfGeom.getSpacing(); // degrees
  const int n2 = static_cast<int>(odfDims[0]);
  const int nPHI = static_cast<int>(odfDims[1]);
  const int n1 = static_cast<int>(odfDims[2]);

  // 2. Reconstruct ODFComponents from the selected Float64 cell arrays. The
  // store is row-major ZYX (phi1 slowest, phi2 fastest), exactly the flat
  // layout gridToODFComponent expects, so we copy values in their natural
  // order.
  std::vector<mtrsim::ODFComponent> components;
  components.reserve(m_InputValues->odf_component_paths.size());
  for(const auto& path : m_InputValues->odf_component_paths)
  {
    const auto& arr = m_DataStructure.getDataRefAs<Float64ArrayType>(path);
    const auto& store = arr.getDataStoreRef();
    std::vector<double> values(store.begin(), store.end());
    components.push_back(mtrsim::GridToOdfComponent(values, n1, nPHI, n2, static_cast<double>(odfSpacing[2]), static_cast<double>(odfSpacing[1]), static_cast<double>(odfSpacing[0])));
  }

  // 3. SimulationParams — always built from m_InputValues, which executeImpl has
  // already populated from either the config file or the manual UI fields.
  mtrsim::SimulationParams params;
  params.x_len = m_InputValues->physical_size[0];
  params.y_len = m_InputValues->physical_size[1];
  params.z_len = m_InputValues->physical_size[2];
  params.dx = m_InputValues->physical_spacing[0];
  params.dy = m_InputValues->physical_spacing[1];
  params.dz = m_InputValues->physical_spacing[2];
  params.volume_fractions = m_InputValues->volume_fractions[0]; // 1 row
  params.theta_list = m_InputValues->theta_list;
  params.seed = m_InputValues->seed;
  // The rng below is seeded from m_InputValues->seed; SimulationParams::seed is
  // not consulted by simulateMTR, but we populate it for completeness.

  if(m_ShouldCancel)
  {
    return {};
  }

  // 4. Run simulation (SIMPLNX z,y,x order out).
  m_MessageHandler.sendInfoMessage("Running MTR simulation (this may take a while for large volumes)...");
  std::mt19937_64 rng(m_InputValues->seed);
  FilterObserver observer{m_MessageHandler, m_ShouldCancel};
  mtrsim::MTRSimResult sim;
  try
  {
    sim = mtrsim::SimulateMtr(params, components, rng, n1, nPHI, n2, &observer);
  } catch(const std::exception& e)
  {
    return MakeErrorResult(-13550, fmt::format("MTR simulation failed: {}", e.what()));
  }

  if(m_ShouldCancel || sim.cancelled)
  {
    return {};
  }

  // 5. Write MTR ids + Euler arrays (output geometry cell AM).
  const DataPath cellAm = m_InputValues->output_geometry_path.createChildPath(m_InputValues->cell_attr_mat_name);
  auto& mtrIds = m_DataStructure.getDataRefAs<Int32ArrayType>(cellAm.createChildPath(m_InputValues->mtr_ids_array_name));
  auto& eulers = m_DataStructure.getDataRefAs<Float32ArrayType>(cellAm.createChildPath(m_InputValues->eulers_array_name));
  auto& mtrStore = mtrIds.getDataStoreRef();
  auto& eulerStore = eulers.getDataStoreRef();

  const std::size_t n = sim.mtr_index.size();
  for(std::size_t i = 0; i < n; ++i)
  {
    mtrStore[i] = sim.mtr_index[i];
    eulerStore[i * 3 + 0] = static_cast<float>(sim.phi1[i]);
    eulerStore[i * 3 + 1] = static_cast<float>(sim.phi[i]);
    eulerStore[i * 3 + 2] = static_cast<float>(sim.phi2[i]);
  }

  m_MessageHandler.sendInfoMessage("MTR simulation complete.");

  if(m_InputValues->generate_polar_coloring)
  {
    m_MessageHandler.sendInfoMessage("Computing polar coloring...");
    auto colorResult = applyPolarColoring(sim, cellAm);
    if(colorResult.invalid())
    {
      return colorResult;
    }
  }

  return {};
}

// -----------------------------------------------------------------------------
Result<> MTRSim::applyPolarColoring(const mtrsim::MTRSimResult& sim, const DataPath& cellAttrMatPath)
{
  const std::size_t n = sim.phi1.size();
  Eigen::VectorXd phi1 = Eigen::Map<const Eigen::VectorXd>(sim.phi1.data(), static_cast<Eigen::Index>(n));
  Eigen::VectorXd phi = Eigen::Map<const Eigen::VectorXd>(sim.phi.data(), static_cast<Eigen::Index>(n));
  Eigen::VectorXd phi2 = Eigen::Map<const Eigen::VectorXd>(sim.phi2.data(), static_cast<Eigen::Index>(n));

  mtrsim::IPFMapper mapper{mtrsim::CrystalSystem::HCP};
  const std::vector<mtrsim::RGBColor> colors = mapper.eulerToColors(phi1, phi, phi2, {0.0, 0.0, 1.0}, mtrsim::IPFColorScheme::MatLab);

  auto& rgb = m_DataStructure.getDataRefAs<UInt8ArrayType>(cellAttrMatPath.createChildPath(m_InputValues->polar_colors_array_name));
  auto& store = rgb.getDataStoreRef();
  for(std::size_t i = 0; i < n; ++i)
  {
    store[i * 3 + 0] = colors[i].r;
    store[i * 3 + 1] = colors[i].g;
    store[i * 3 + 2] = colors[i].b;
  }
  return {};
}

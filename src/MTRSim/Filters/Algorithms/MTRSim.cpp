#include "MTRSim.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"

#include "LibMTRSim/MTRSimDriver.hpp"

#include <fmt/format.h>

#include <cmath>
#include <exception>
#include <random>
#include <vector>

using namespace nx::core;

// -----------------------------------------------------------------------------
MTRSim::MTRSim(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, MTRSimInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

// -----------------------------------------------------------------------------
MTRSim::~MTRSim() noexcept = default;

// -----------------------------------------------------------------------------
Result<> MTRSim::operator()()
{
  m_MessageHandler(IFilter::Message::Type::Info, "Reading ODF geometry...");

  // 1. ODF grid geometry (degrees) from the input ODF ImageGeom.
  const auto& odfGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->inputOdfGeometryPath);
  const SizeVec3 odfDims = odfGeom.getDimensions();  // X=phi2, Y=PHI, Z=phi1
  const FloatVec3 odfSpacing = odfGeom.getSpacing(); // degrees
  const int n2 = static_cast<int>(odfDims[0]);
  const int nPHI = static_cast<int>(odfDims[1]);
  const int n1 = static_cast<int>(odfDims[2]);

  // 2. Reconstruct ODFComponents from the selected Float64 cell arrays. The
  // store is row-major ZYX (phi1 slowest, phi2 fastest), exactly the flat
  // layout gridToODFComponent expects, so we copy values in their natural order.
  std::vector<mtrsim::ODFComponent> components;
  components.reserve(m_InputValues->odfComponentPaths.size());
  for(const auto& path : m_InputValues->odfComponentPaths)
  {
    const auto& arr = m_DataStructure.getDataRefAs<Float64Array>(path);
    const auto& store = arr.getDataStoreRef();
    std::vector<double> values(store.begin(), store.end());
    components.push_back(mtrsim::gridToODFComponent(values, n1, nPHI, n2, static_cast<double>(odfSpacing[2]), static_cast<double>(odfSpacing[1]), static_cast<double>(odfSpacing[0])));
  }

  // 3. SimulationParams from filter inputs.
  mtrsim::SimulationParams params;
  params.xLen = m_InputValues->physicalSize[0];
  params.yLen = m_InputValues->physicalSize[1];
  params.zLen = m_InputValues->physicalSize[2];
  params.dx = m_InputValues->physicalSpacing[0];
  params.dy = m_InputValues->physicalSpacing[1];
  params.dz = m_InputValues->physicalSpacing[2];
  params.volumeFractions = m_InputValues->volumeFractions[0]; // 1 row
  params.thetaList = m_InputValues->thetaList;
  params.seed = m_InputValues->seed;

  if(m_ShouldCancel)
  {
    return {};
  }

  // 4. Run simulation (SIMPLNX z,y,x order out).
  m_MessageHandler(IFilter::Message::Type::Info, "Running MTR simulation (this may take a while for large volumes)...");
  std::mt19937_64 rng(m_InputValues->seed);
  mtrsim::MTRSimResult sim;
  try
  {
    sim = mtrsim::simulateMTR(params, components, rng, n1, nPHI, n2);
  } catch(const std::exception& e)
  {
    return MakeErrorResult(-13050, fmt::format("MTR simulation failed: {}", e.what()));
  }

  if(m_ShouldCancel)
  {
    return {};
  }

  // 5. Write MTR ids + Euler arrays (output geometry cell AM).
  const DataPath cellAm = m_InputValues->outputGeometryPath.createChildPath(m_InputValues->cellAttrMatName);
  auto& mtrIds = m_DataStructure.getDataRefAs<Int32Array>(cellAm.createChildPath(m_InputValues->mtrIdsArrayName));
  auto& eulers = m_DataStructure.getDataRefAs<Float32Array>(cellAm.createChildPath(m_InputValues->eulersArrayName));
  auto& mtrStore = mtrIds.getDataStoreRef();
  auto& eulerStore = eulers.getDataStoreRef();

  const std::size_t N = sim.mtrIndex.size();
  for(std::size_t i = 0; i < N; ++i)
  {
    mtrStore[i] = sim.mtrIndex[i];
    eulerStore[i * 3 + 0] = static_cast<float>(sim.phi1[i]);
    eulerStore[i * 3 + 1] = static_cast<float>(sim.phi[i]);
    eulerStore[i * 3 + 2] = static_cast<float>(sim.phi2[i]);
  }

  m_MessageHandler(IFilter::Message::Type::Info, "MTR simulation complete.");
  // Polar coloring (Task 8) is handled in a later task; the optional array (if
  // requested) is left created-but-unfilled here and populated next task.
  return {};
}

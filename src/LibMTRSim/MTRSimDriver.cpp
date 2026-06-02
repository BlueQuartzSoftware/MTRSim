#include "MTRSimDriver.hpp"

#include "ODFSampler.hpp"
#include "PGRFSimulation.hpp"

#include <algorithm>
#include <cmath>
#include <fmt/format.h>
#include <numbers>
#include <stdexcept>

namespace mtrsim
{

ODFComponent buildUniformODF(int n1, int nPHI, int n2)
{
  const int nTotal = n1 * nPHI * n2;
  const double twoPiOverN1 = 2.0 * std::numbers::pi / static_cast<double>(n1);
  const double piOverNPHI = std::numbers::pi / static_cast<double>(nPHI);
  const double twoPiOverN2 = 2.0 * std::numbers::pi / static_cast<double>(n2);

  Eigen::VectorXd phi1Bins(nTotal);
  Eigen::VectorXd phiBins(nTotal);
  Eigen::VectorXd phi2Bins(nTotal);

  for(int ix = 0; ix < nTotal; ++ix)
  {
    const int i1 = ix / (nPHI * n2);
    const int iPHI = (ix % (nPHI * n2)) / n2;
    const int i2 = ix % n2;
    phi1Bins[ix] = (i1 + 0.5) * twoPiOverN1;
    phiBins[ix] = (iPHI + 0.5) * piOverNPHI;
    phi2Bins[ix] = (i2 + 0.5) * twoPiOverN2;
  }

  ODFComponent uni;
  uni.odfVal = Eigen::VectorXd::Constant(nTotal, 1.0 / static_cast<double>(nTotal));
  uni.phi1Bins = std::move(phi1Bins);
  uni.phiBins = std::move(phiBins);
  uni.phi2Bins = std::move(phi2Bins);
  return uni;
}

ODFComponent gridToODFComponent(const std::vector<double>& values, int n1, int nPHI, int n2, double stepDeg1, double stepDegPHI, double stepDeg2)
{
  const int nTotal = n1 * nPHI * n2;
  const double deg2rad = std::numbers::pi / 180.0;
  const double s1 = stepDeg1 * deg2rad;
  const double sP = stepDegPHI * deg2rad;
  const double s2 = stepDeg2 * deg2rad;

  Eigen::VectorXd phi1Bins(nTotal);
  Eigen::VectorXd phiBins(nTotal);
  Eigen::VectorXd phi2Bins(nTotal);
  for(int ix = 0; ix < nTotal; ++ix)
  {
    const int i1 = ix / (nPHI * n2);
    const int iPHI = (ix % (nPHI * n2)) / n2;
    const int i2 = ix % n2;
    phi1Bins[ix] = (i1 + 0.5) * s1;
    phiBins[ix] = (iPHI + 0.5) * sP;
    phi2Bins[ix] = (i2 + 0.5) * s2;
  }

  ODFComponent c;
  c.odfVal = Eigen::Map<const Eigen::VectorXd>(values.data(), nTotal);
  const double total = c.odfVal.sum();
  if(total > 0.0)
  {
    c.odfVal /= total;
  }
  c.phi1Bins = std::move(phi1Bins);
  c.phiBins = std::move(phiBins);
  c.phi2Bins = std::move(phi2Bins);
  return c;
}

MTRSimResult simulateMTR(const SimulationParams& params, const std::vector<ODFComponent>& odfComponents, std::mt19937_64& rng, int n1, int nPHI, int n2, ISimulationObserver* observer)
{
  const int nx = static_cast<int>(std::round(params.xLen / params.dx));
  const int ny = static_cast<int>(std::round(params.yLen / params.dy));
  const int nz = std::max(static_cast<int>(std::round(params.zLen / params.dz)), 1);
  const int N = nx * ny * nz;

  if(static_cast<int>(odfComponents.size()) != static_cast<int>(params.volumeFractions.size()))
  {
    throw std::invalid_argument("simulateMTR: odfComponents count must equal volumeFractions count");
  }

  auto cancelled = [&]() { return observer != nullptr && observer->shouldCancel(); };
  auto report = [&](int64_t done, int64_t total, const std::string& msg) {
    if(observer != nullptr)
    {
      observer->updateProgress(done, total, msg);
    }
  };

  // 1. PGRF assignment (sim-ordered, 1-based component ids).
  report(0, 100, "Running plurigaussian field simulation");
  PGRFSimulation pgrf{rng};
  const PGRFResult pgrf_result = pgrf.run(params, observer); // throws on bad dims
  if(cancelled())
  {
    MTRSimResult out;
    out.cancelled = true;
    return out;
  }

  if(static_cast<int>(pgrf_result.mtrIndex.size()) != N)
  {
    throw std::runtime_error("simulateMTR: PGRF result size does not match grid dimensions");
  }

  // 2. Sample N orientations per component against the uniform reference.
  const ODFComponent uniformOdf = buildUniformODF(n1, nPHI, n2);
  const int numComponents = static_cast<int>(odfComponents.size());
  std::vector<Eigen::MatrixXd> orientSamples(static_cast<std::size_t>(numComponents));
  ODFSampler sampler{rng};
  for(int j = 0; j < numComponents; ++j)
  {
    report(j, numComponents, fmt::format("Sampling orientations (component {}/{})", j + 1, numComponents));
    orientSamples[static_cast<std::size_t>(j)] = sampler.sampleN(N, odfComponents[static_cast<std::size_t>(j)], uniformOdf, observer);
    if(cancelled())
    {
      MTRSimResult out;
      out.cancelled = true;
      return out;
    }
  }

  // 3. Assign per-voxel orientation by component (sim order).
  std::vector<int32_t> mtrSim(static_cast<std::size_t>(N));
  std::vector<double> phi1Sim(static_cast<std::size_t>(N));
  std::vector<double> phiSim(static_cast<std::size_t>(N));
  std::vector<double> phi2Sim(static_cast<std::size_t>(N));
  for(int i = 0; i < N; ++i)
  {
    const int comp = pgrf_result.mtrIndex[i] - 1;
    mtrSim[static_cast<std::size_t>(i)] = pgrf_result.mtrIndex[i];
    phi1Sim[static_cast<std::size_t>(i)] = orientSamples[static_cast<std::size_t>(comp)](i, 0);
    phiSim[static_cast<std::size_t>(i)] = orientSamples[static_cast<std::size_t>(comp)](i, 1);
    phi2Sim[static_cast<std::size_t>(i)] = orientSamples[static_cast<std::size_t>(comp)](i, 2);
  }

  // 4. Remap to SIMPLNX z,y,x order.
  report(100, 100, "Finalizing microstructure");
  MTRSimResult out;
  out.nx = nx;
  out.ny = ny;
  out.nz = nz;
  out.mtrIndex = remapSimToZYX(mtrSim, nx, ny, nz);
  out.phi1 = remapSimToZYX(phi1Sim, nx, ny, nz);
  out.phi = remapSimToZYX(phiSim, nx, ny, nz);
  out.phi2 = remapSimToZYX(phi2Sim, nx, ny, nz);
  return out;
}

} // namespace mtrsim

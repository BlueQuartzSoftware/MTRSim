// Port of PGRF_simulation.m — plurigaussian random field simulation
// orchestrator. Daniel M. Sparkman, 08/28/2013.

#include "PGRFSimulation.hpp"

#include "AssignmentRule.hpp"
#include "GPGenerator.hpp"

#include <algorithm>
#include <cmath>
#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <stdexcept>

namespace mtrsim
{

PGRFSimulation::PGRFSimulation(std::mt19937_64& rng)
: m_Rng(rng)
{
}

// ─────────────────────────────────────────────────────────────────────────────
// run
//
// Orchestrates the full plurigaussian simulation:
//   1. Derive grid dimensions from params.
//   2. Select Gaussian thresholds (AssignmentRule::selectThresholds).
//   3. For each latent Gaussian field, call GPGenerator::generate.
//   4. Classify every voxel (AssignmentRule::evaluate).
//
// The correlation model is exponential and anisotropic (non-periodic), matching
// the hardcoded MATLAB settings:
//   corr_func_name            = 'exp'
//   correlation_function_selected = 'anisotropic'
//   boundary_conditions       = 'nonperiodic'
//   mean_function_selected    = 'stationary'  (mu = 0 for all fields)

PGRFResult PGRFSimulation::run(const SimulationParams& params, ISimulationObserver* observer)
{
  // ── Grid dimensions ──────────────────────────────────────────────────────
  const int nx = static_cast<int>(std::round(params.xLen / params.dx));
  const int ny = static_cast<int>(std::round(params.yLen / params.dy));
  const int nz = std::max(static_cast<int>(std::round(params.zLen / params.dz)), 1);
  const int N = nx * ny * nz;
  const int numComponents = static_cast<int>(params.volumeFractions.size());
  const int numGaussians = numComponents - 1;

  if(numGaussians < 1)
  {
    throw std::invalid_argument("PGRFSimulation: need at least 2 volume fraction components");
  }
  if(static_cast<int>(params.thetaList.size()) < numGaussians)
  {
    throw std::invalid_argument("PGRFSimulation: thetaList must have one row per latent Gaussian");
  }

  spdlog::info("PGRFSimulation: grid {}x{}x{} = {} voxels, {} components, {} "
               "latent fields",
               nx, ny, nz, N, numComponents, numGaussians);

  // ── Select assignment-rule thresholds ────────────────────────────────────
  spdlog::info("PGRFSimulation: selecting assignment rule thresholds...");
  AssignmentRule ar(m_Rng);
  const AssignmentRuleThresholds thresholds = ar.selectThresholds(params.volumeFractions);

  // ── Exponential correlation function: rho(lag, theta) = exp(-|lag|/theta) ─
  // Matches MATLAB: corr_func_name = 'exp', boundary_conditions = 'nonperiodic'
  auto corrFn = [](double lag, double theta) -> double { return std::exp(-std::abs(lag) / theta); };

  GPGenerator gpGen(m_Rng, corrFn);

  // ── Simulate each latent Gaussian field ──────────────────────────────────
  // Stationary mean: mu_const = 0 for all fields (matches MATLAB default)
  Eigen::MatrixXd zAll = Eigen::MatrixXd::Zero(N, numGaussians);

  for(int h = 0; h < numGaussians; ++h)
  {
    if(observer != nullptr)
    {
      if(observer->shouldCancel())
      {
        return PGRFResult{}; // empty; caller checks observer->shouldCancel()
      }
      observer->updateProgress(h, numGaussians, fmt::format("Simulating latent Gaussian field {}/{}", h + 1, numGaussians));
    }
    spdlog::info("PGRFSimulation: simulating latent Gaussian Y{} ...", h + 1);

    const auto& thetaRow = params.thetaList[static_cast<std::size_t>(h)];
    if(thetaRow.size() < 3)
    {
      throw std::invalid_argument("PGRFSimulation: each thetaList row must have 3 elements [theta_x, "
                                  "theta_y, theta_z]");
    }

    const std::array<double, 3> theta = {thetaRow[0], thetaRow[1], thetaRow[2]};
    zAll.col(h) = gpGen.generate(params.dx, params.dy, params.dz, theta, nx, ny, nz);
  }

  // ── Apply assignment rule ─────────────────────────────────────────────────
  spdlog::info("PGRFSimulation: applying assignment rule...");
  const Eigen::VectorXi mtrIndex = ar.evaluate(zAll, thresholds);

  // Log empirical volume fractions for verification
  for(int j = 1; j <= numComponents; ++j)
  {
    const int count = (mtrIndex.array() == j).count();
    spdlog::info("  P{} empirical = {:.3f}  (target {:.3f})", j, static_cast<double>(count) / static_cast<double>(N), params.volumeFractions[static_cast<std::size_t>(j - 1)]);
  }

  return PGRFResult{mtrIndex, zAll};
}

} // namespace mtrsim

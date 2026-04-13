#pragma once

#include "SimulationParams.hpp"
#include "mtrsim_export.h"
#include <Eigen/Dense>
#include <random>

namespace mtrsim
{

/**
 * @brief Result of the plurigaussian random field simulation.
 */
struct MTRSIM_EXPORT PGRFResult
{
  // Component assignment for each voxel (1-based index), length N
  Eigen::VectorXi mtrIndex;

  // Latent Gaussian field values, shape [N x numGaussians]
  Eigen::MatrixXd latentFields;
};

/**
 * @brief Simulates the plurigaussian random field (PGRF) that assigns each
 *        voxel to an MTR component class.
 *
 * Equivalent to PGRF_simulation.m.  Orchestrates GPGenerator and
 * AssignmentRule to produce a categorical assignment over the simulation
 * volume.
 */
class MTRSIM_EXPORT PGRFSimulation
{
public:
  explicit PGRFSimulation(std::mt19937_64& rng);

  /**
   * @brief Run the PGRF simulation.
   *
   * @param params  Fully populated SimulationParams
   * @return        PGRFResult containing voxel assignments and latent fields
   */
  PGRFResult run(const SimulationParams& params);

private:
  std::mt19937_64& m_Rng;
};

} // namespace mtrsim

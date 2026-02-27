#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace mtrsim
{

/**
 * @brief Holds all parameters that drive a single MTR simulation run.
 *
 * These map directly to the top-level parameter block in simulate_MTRs.m.
 */
struct SimulationParams
{
  // Volume dimensions [mm]
  double xLen = 1.5 * 25.4;
  double yLen = 0.5 * 25.4;
  double zLen = 0.0;

  // Voxel spacing [mm]
  double dx = 0.02;
  double dy = 0.02;
  double dz = 0.02;

  // Volume fractions for each MTR component ODF (must sum to 1)
  std::vector<double> volumeFractions = {0.30, 0.35, 0.35};

  // Correlation length parameters [num_gaussians x 3 (x,y,z directions)]
  // Stored row-major: theta_list(gaussian, direction)
  std::vector<std::vector<double>> thetaList = {
      {0.10, 0.45, 0.10},
      {0.08, 0.37, 0.08},
  };

  // Nugget variance per component (currently unused in simulation)
  std::vector<double> nuggetVariance = {0.67, 0.71, 0.72};

  // Paths
  std::string odfInputPath = "data/simulation_ODF.h5";
  std::string outputDir = ".";

  // Random seed (0 = use std::random_device)
  uint64_t seed = 0;
};

} // namespace mtrsim

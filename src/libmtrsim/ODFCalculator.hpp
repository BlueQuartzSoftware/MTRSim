#pragma once

#include "ODFSampler.hpp"
#include "mtrsim_export.h"
#include <Eigen/Dense>

namespace mtrsim {

/**
 * @brief Computes a discrete ODF histogram from a set of Euler angles.
 *
 * Equivalent to calc_ODF.m.  Applies crystal symmetry operations (HCP by
 * default) and optional nearest-neighbour smoothing before binning.
 */
class MTRSIM_EXPORT ODFCalculator {
public:
  ODFCalculator() = default;

  /**
   * @brief Compute the ODF from a collection of orientations.
   *
   * @param phi1     phi1 angles [rad], length N
   * @param phi      PHI  angles [rad], length N
   * @param phi2     phi2 angles [rad], length N
   * @param degSpacing  Bin width in degrees (default 5°)
   * @return         Populated ODFComponent with normalised ODFval
   */
  ODFComponent compute(const Eigen::VectorXd &phi1, const Eigen::VectorXd &phi,
                       const Eigen::VectorXd &phi2, double degSpacing = 5.0);
};

} // namespace mtrsim

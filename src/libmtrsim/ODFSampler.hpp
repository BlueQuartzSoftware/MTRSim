#pragma once

#include "mtrsim_export.h"
#include <Eigen/Dense>
#include <random>
#include <vector>

namespace mtrsim {

/**
 * @brief Discrete ODF data for one component.
 *
 * Matches the MATLAB struct fields stored in simulation_ODF.h5.
 */
struct MTRSIM_EXPORT ODFComponent {
  Eigen::VectorXd
      odfVal; ///< Probability mass for each Euler-space bin (N_bins,)
  Eigen::VectorXd phi1Bins; ///< phi1 bin centres [rad]
  Eigen::VectorXd phiBins;  ///< PHI  bin centres [rad]
  Eigen::VectorXd phi2Bins; ///< phi2 bin centres [rad]
};

/**
 * @brief Sampled orientation (Bunge Euler angles, radians).
 */
struct MTRSIM_EXPORT EulerAngles {
  double phi1 = 0.0;
  double phi = 0.0;
  double phi2 = 0.0;
};

/**
 * @brief Draws orientations from a discrete ODF by inverse-CDF sampling.
 *
 * Combines sample_orientation_from_ODF.m and
 * sample_N_orientations_from_ODF.m.
 */
class MTRSIM_EXPORT ODFSampler {
public:
  explicit ODFSampler(std::mt19937_64 &rng);

  /**
   * @brief Draw N orientations from the given component ODF.
   *
   * @param n         Number of orientations to draw
   * @param component Target ODF component
   * @param uniform   Uniform (reference) ODF component for bin coordinates
   * @return          Matrix of shape [N x 3]: columns are phi1, PHI, phi2 [rad]
   */
  Eigen::MatrixXd sampleN(int n, const ODFComponent &component,
                          const ODFComponent &uniform);

  /**
   * @brief Draw a single orientation from the given component ODF.
   */
  EulerAngles sampleOne(const ODFComponent &component,
                        const ODFComponent &uniform);

private:
  std::mt19937_64 &m_Rng;
};

} // namespace mtrsim

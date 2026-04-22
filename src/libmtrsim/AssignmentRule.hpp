#pragma once

#include "mtrsim_export.h"
#include <Eigen/Dense>
#include <random>
#include <vector>

namespace mtrsim {

/**
 * @brief Result of threshold selection for the plurigaussian assignment rule.
 */
struct MTRSIM_EXPORT AssignmentRuleThresholds {
  int numGaussians = 0;

  // min_thresholds(component, gaussian)  — shape [numComponents x numGaussians]
  Eigen::MatrixXd minThresholds;

  // max_thresholds(component, gaussian)  — shape [numComponents x numGaussians]
  Eigen::MatrixXd maxThresholds;
};

/**
 * @brief Selects and evaluates the plurigaussian assignment rule.
 *
 * Combines the logic of select_AR.m and eval_AR.m.
 */
class MTRSIM_EXPORT AssignmentRule {
public:
  explicit AssignmentRule(std::mt19937_64 &rng);

  /**
   * @brief Determine Gaussian thresholds that yield the requested volume
   *        fractions P.
   *
   * @param volumeFractions  Target probability for each MTR component
   * @return Threshold matrices for use in evaluate()
   */
  AssignmentRuleThresholds
  selectThresholds(const std::vector<double> &volumeFractions);

  /**
   * @brief Classify each voxel given realizations of the latent Gaussians.
   *
   * @param z           Matrix of latent field values, shape [N x numGaussians]
   * @param thresholds  Thresholds returned by selectThresholds()
   * @return            Integer component index (1-based) for each voxel, length
   * N
   */
  Eigen::VectorXi evaluate(const Eigen::MatrixXd &z,
                           const AssignmentRuleThresholds &thresholds) const;

private:
  std::mt19937_64 &m_Rng;
};

} // namespace mtrsim

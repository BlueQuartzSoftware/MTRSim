#pragma once

#include <Eigen/Dense>
#include <random>
#include <utility>

namespace mtrsim
{

/**
 * @brief Quasi-Monte Carlo estimator for the multivariate normal CDF.
 *
 * Port of qsimvn.m (Alan Genz, 1992).  Estimates
 *
 *   P = Prob( a <= X <= b )
 *
 * where X ~ N(0, R) using a randomised quasi-random rule with m sample
 * points.
 *
 * License note: the underlying algorithm is due to Alan Genz (BSD-style).
 * Full attribution is preserved in matlab/qsimvn.m.
 */
class QSimVN
{
public:
  /**
   * @param rng Seeded RNG engine (shared with the rest of the simulation).
   */
  explicit QSimVN(std::mt19937_64& rng);

  /**
   * @brief Estimate the MVN probability.
   *
   * @param m   Number of quasi-random sample points
   * @param r   Positive-definite covariance matrix (n×n)
   * @param a   Lower integration limits (length n, may be -inf)
   * @param b   Upper integration limits (length n, may be +inf)
   * @return    {probability estimate, error estimate}
   */
  std::pair<double, double> compute(int m,
                                    const Eigen::MatrixXd& r,
                                    const Eigen::VectorXd& a,
                                    const Eigen::VectorXd& b);

private:
  // Cholesky decomposition with reordering (chlrdr in original)
  struct ChlrdrResult
  {
    Eigen::MatrixXd ch;
    Eigen::VectorXd ap;
    Eigen::VectorXd bp;
  };

  ChlrdrResult chlrdr(const Eigen::MatrixXd& r,
                      const Eigen::VectorXd& a,
                      const Eigen::VectorXd& b) const;

  double mvndns(int n, const Eigen::MatrixXd& ch,
                double ci, double dci,
                const Eigen::VectorXd& x,
                const Eigen::VectorXd& a,
                const Eigen::VectorXd& b) const;

  static double phi(double z);
  static double phiInv(double p);

  std::mt19937_64& m_Rng;
};

} // namespace mtrsim

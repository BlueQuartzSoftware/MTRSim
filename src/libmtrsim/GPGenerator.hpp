#pragma once

#include "mtrsim_export.h"
#include <Eigen/Dense>
#include <functional>
#include <random>

namespace mtrsim {

/**
 * @brief Generates a realization of a separable Gaussian random field via
 *        the Kronecker-Cholesky method.
 *
 * Equivalent to gp_generator.m.
 *
 * The covariance in each spatial direction is built from a user-supplied
 * correlation function rho(lag, theta).  The field is drawn using the property:
 *
 *   kron(Rx, kron(Ry, Rz)) * vec(Z) = Rz * reshape(Z, nx*ny, nz)' * ...
 *
 * where Rx, Ry, Rz are the Cholesky factors of the per-direction covariance
 * matrices.
 */
class MTRSIM_EXPORT GPGenerator {
public:
  using CorrelationFn = std::function<double(double lag, double theta)>;

  /**
   * @param rng    Seeded Mersenne-Twister engine (passed by reference so the
   *               caller controls the global RNG state).
   * @param corrFn Correlation function: rho(|lag|, theta) → [0, 1]
   */
  explicit GPGenerator(std::mt19937_64 &rng, CorrelationFn corrFn);

  /**
   * @brief Draw one realization of the GP on an nx × ny × nz grid.
   *
   * @param hx     Voxel spacing in x
   * @param hy     Voxel spacing in y
   * @param hz     Voxel spacing in z
   * @param theta  Correlation length parameters [theta_x, theta_y, theta_z]
   * @param nx     Number of voxels in x
   * @param ny     Number of voxels in y
   * @param nz     Number of voxels in z
   * @return       Flattened field of length nx*ny*nz (z-major ordering)
   */
  Eigen::VectorXd generate(double hx, double hy, double hz,
                           const std::array<double, 3> &theta, int nx, int ny,
                           int nz);

private:
  Eigen::MatrixXd buildCovarianceMatrix(int n, double spacing,
                                        double theta) const;

  std::mt19937_64 &m_Rng;
  CorrelationFn m_CorrFn;
};

} // namespace mtrsim

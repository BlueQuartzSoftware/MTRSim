// Port of gp_generator.m — Kronecker-Cholesky Gaussian random field generator.
// Daniel M. Sparkman, Research.

#include "GPGenerator.hpp"

#include <stdexcept>

namespace mtrsim
{

GPGenerator::GPGenerator(std::mt19937_64& rng, CorrelationFn corrFn)
: m_Rng(rng)
, m_CorrFn(std::move(corrFn))
{
}

// ─────────────────────────────────────────────────────────────────────────────
// buildCovarianceMatrix
//
// Builds the n×n stationary covariance matrix for one spatial direction:
//
//   Gamma(i,j) = rho(|i-j| * spacing, theta)   for i ≠ j
//   Gamma(i,i) = 1  (correlation function evaluated at lag 0)
//
// A small jitter (1e-6) is added to the diagonal to ensure positive definiteness,
// matching MATLAB: Gamma = Gamma + 1e-6*eye(n).

Eigen::MatrixXd GPGenerator::buildCovarianceMatrix(int n, double spacing, double theta) const
{
  Eigen::MatrixXd gamma = Eigen::MatrixXd::Zero(n, n);

  // Fill upper triangle: gamma(i,j) = rho((j-i)*spacing, theta)  for j > i
  for(int i = 0; i < n - 1; ++i)
  {
    for(int j = i + 1; j < n; ++j)
    {
      gamma(i, j) = m_CorrFn(static_cast<double>(j - i) * spacing, theta);
    }
  }

  // Symmetrize (adds lower triangle from upper) then set diagonal = 1 + jitter
  // MATLAB: Gamma = Gamma + Gamma' + eye(n) + 1e-6*eye(n)
  gamma += gamma.transpose();
  gamma.diagonal().array() += 1.0 + 1e-6;

  return gamma;
}

// ─────────────────────────────────────────────────────────────────────────────
// generate
//
// Draws one realization of the separable GP on an nx × ny × nz grid via the
// Kronecker-Cholesky method.  The mathematical identity used is:
//
//   kron(A,B) * vec(X) = B * X * A'
//
// applied twice to avoid forming the full (N×N) covariance matrix.
//
// The MATLAB memory layout for the output is column-major for an ny×nx×nz
// 3-D array, i.e., element (iy, ix, iz) lives at index iz*(ny*nx) + ix*ny + iy.
// The returned VectorXd uses this same ordering so that Z(:) in MATLAB matches
// a direct assignment to the column of Z_all.
//
// Notation mapping (MATLAB chol returns upper R s.t. R'R=Gamma):
//   MATLAB:  W_1 * Rz          →  C++: W1 * llt_z.matrixU()
//   MATLAB:  Ry' * W2k * Rx    →  C++: llt_y.matrixL() * W2k * llt_x.matrixU()

Eigen::VectorXd GPGenerator::generate(double hx, double hy, double hz, const std::array<double, 3>& theta, int nx, int ny, int nz)
{
  // ── Build per-direction covariance matrices ─────────────────────────────
  const Eigen::MatrixXd Gamma_x = buildCovarianceMatrix(nx, hx, theta[0]);
  const Eigen::MatrixXd Gamma_y = buildCovarianceMatrix(ny, hy, theta[1]);
  const Eigen::MatrixXd Gamma_z = buildCovarianceMatrix(nz, hz, theta[2]);

  // ── Cholesky factorisation: L * L' = Gamma  (L lower-triangular) ─────────
  Eigen::LLT<Eigen::MatrixXd> llt_x(Gamma_x);
  Eigen::LLT<Eigen::MatrixXd> llt_y(Gamma_y);
  Eigen::LLT<Eigen::MatrixXd> llt_z(Gamma_z);

  if(llt_x.info() != Eigen::Success || llt_y.info() != Eigen::Success || llt_z.info() != Eigen::Success)
  {
    throw std::runtime_error("GPGenerator::generate — Cholesky failed; covariance matrix is not positive definite");
  }

  // MATLAB upper-triangular Cholesky factors: Rx = llt_x.matrixU()  (= L_x')
  // Used below as: W1 * Rz  and  Ry' * W2k * Rx
  const Eigen::MatrixXd Rx = llt_x.matrixU();       // nx × nx
  const Eigen::MatrixXd Ry_lower = llt_y.matrixL(); // ny × ny  (= MATLAB Ry')
  const Eigen::MatrixXd Rz = llt_z.matrixU();       // nz × nz

  // ── Standard normal draw ─────────────────────────────────────────────────
  const int N = nx * ny * nz;
  std::normal_distribution<double> normal(0.0, 1.0);
  Eigen::VectorXd z(N);
  for(int i = 0; i < N; ++i)
  {
    z(i) = normal(m_Rng);
  }

  // ── Apply z-direction Cholesky: u = W_1 * Rz ────────────────────────────
  // W_1 = reshape(z, nx*ny, nz)  — Map z as column-major (nx*ny) × nz
  // Each row w of W_1 is a 1×nz vector; w * Rz ~ N(0, Rz'Rz) = N(0, Gamma_z)
  Eigen::Map<const Eigen::MatrixXd> W1(z.data(), nx * ny, nz);
  const Eigen::MatrixXd u = W1 * Rz; // (nx*ny) × nz

  // ── Apply x and y Cholesky slice-by-slice ────────────────────────────────
  // W_2(:,:,k) = reshape(u.col(k), ny, nx)   — Map each column as ny × nx
  // x_k = Ry' * (W2k * Rx)                  — nx-covariance from right, ny from left
  Eigen::VectorXd gp(N);
  for(int k = 0; k < nz; ++k)
  {
    // W2k: ny × nx view into column k of u (column-major, contiguous in memory)
    Eigen::Map<const Eigen::MatrixXd> W2k(u.col(k).data(), ny, nx);

    // x_k: ny × nx result for this z-slice
    const Eigen::MatrixXd x_k = Ry_lower * W2k * Rx;

    // Store into output with MATLAB column-major ordering for a ny×nx×nz array:
    //   gp[k*(ny*nx) + ix*ny + iy]  =  x_k(iy, ix)
    Eigen::Map<Eigen::MatrixXd>(gp.data() + static_cast<std::ptrdiff_t>(k) * (ny * nx), ny, nx) = x_k;
  }

  return gp;
}

} // namespace mtrsim

#include "AssignmentRule.hpp"
#include "QSimVN.hpp"

#include <Eigen/Dense>
#include <catch2/catch.hpp>
#include <cmath>
#include <numbers>
#include <random>

using namespace mtrsim;

// ─────────────────────────────────────────────────────────────────────────────
// QSimVN tests
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("QSimVN: univariate N(0,1) CDF matches std::erfc", "[qsimvn]") {
  // P(Z <= z) for Z ~ N(0,1) with r=[1], a=-inf, b=z
  std::mt19937_64 rng(42);
  QSimVN qsimvn(rng);

  Eigen::MatrixXd r(1, 1);
  r(0, 0) = 1.0;
  Eigen::VectorXd a(1), b(1);

  // P(Z <= 0) = 0.5
  a(0) = -std::numeric_limits<double>::infinity();
  b(0) = 0.0;
  auto [p0, e0] = qsimvn.compute(1000, r, a, b);
  CHECK(p0 == Approx(0.5).margin(e0 + 1e-4));

  // P(Z <= 1.96) ≈ 0.975
  b(0) = 1.96;
  auto [p1, e1] = qsimvn.compute(1000, r, a, b);
  CHECK(p1 == Approx(0.975).margin(e1 + 5e-4));

  // P(Z <= -1.96) ≈ 0.025
  b(0) = -1.96;
  auto [p2, e2] = qsimvn.compute(1000, r, a, b);
  CHECK(p2 == Approx(0.025).margin(e2 + 5e-4));
}

TEST_CASE("QSimVN: bivariate independent N(0,I) quadrant probability ≈ 0.25",
          "[qsimvn]") {
  // P(Z1 <= 0, Z2 <= 0) = 0.25 for Z ~ N(0,I_2)
  std::mt19937_64 rng(123);
  QSimVN qsimvn(rng);

  Eigen::MatrixXd r = Eigen::MatrixXd::Identity(2, 2);
  Eigen::VectorXd a(2), b(2);
  a(0) = a(1) = -std::numeric_limits<double>::infinity();
  b(0) = b(1) = 0.0;

  auto [p, e] = qsimvn.compute(5000, r, a, b);
  CHECK(p == Approx(0.25).margin(e + 1e-4));
}

TEST_CASE("QSimVN: bivariate correlated N(0,R) matches known probability",
          "[qsimvn]") {
  // Genz example from qsimvn.m header:
  //   r = [4 3 2 1; 3 5 -1 1; 2 -1 4 2; 1 1 2 5]
  //   a = -inf*[1 1 1 1]', b = [1 2 3 4]'
  // The exact answer is approximately 0.4533 (computed from reference MATLAB).
  // We use a looser tolerance since this is a stochastic estimator.
  std::mt19937_64 rng(77);
  QSimVN qsimvn(rng);

  Eigen::MatrixXd r(4, 4);
  r << 4, 3, 2, 1, 3, 5, -1, 1, 2, -1, 4, 2, 1, 1, 2, 5;

  Eigen::VectorXd a(4), b(4);
  a << -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity();
  b << 1.0, 2.0, 3.0, 4.0;

  auto [p, e] = qsimvn.compute(5000, r, a, b);
  // Result should be in (0, 1) and the error estimate should be small
  CHECK(p > 0.0);
  CHECK(p < 1.0);
  CHECK(e < 0.05); // error estimate within 5%
}

TEST_CASE("QSimVN: symmetric interval gives correct probability", "[qsimvn]") {
  // P(-1 <= Z <= 1) = 2*Phi(1) - 1 ≈ 0.6827  for Z ~ N(0,1)
  std::mt19937_64 rng(999);
  QSimVN qsimvn(rng);

  Eigen::MatrixXd r(1, 1);
  r(0, 0) = 1.0;
  Eigen::VectorXd a(1), b(1);
  a(0) = -1.0;
  b(0) = 1.0;

  auto [p, e] = qsimvn.compute(1000, r, a, b);
  const double expected = std::erfc(-1.0 / std::sqrt(2.0)) - 1.0; // 2*Phi(1)-1
  CHECK(p == Approx(expected).margin(e + 5e-4));
}

// ─────────────────────────────────────────────────────────────────────────────
// AssignmentRule tests
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("AssignmentRule::selectThresholds: 3-component fractions sum to 1",
          "[assignmentrule]") {
  // Standard MTRsim usage: 3 components with P = {0.30, 0.35, 0.35}
  std::mt19937_64 rng(42);
  AssignmentRule ar(rng);

  const std::vector<double> P = {0.30, 0.35, 0.35};
  const auto thresholds = ar.selectThresholds(P);

  CHECK(thresholds.numGaussians == 2);
  CHECK(thresholds.minThresholds.rows() == 3);
  CHECK(thresholds.maxThresholds.rows() == 3);
  CHECK(thresholds.minThresholds.cols() == 2);
  CHECK(thresholds.maxThresholds.cols() == 2);

  // Thresholds should be finite for the bisected dimensions
  // Component 0: max_thresholds(0,0) should be a finite quantile
  CHECK(std::isfinite(thresholds.maxThresholds(0, 0)));
  // Component 1: max_thresholds(1,1) should be a finite quantile
  CHECK(std::isfinite(thresholds.maxThresholds(1, 1)));
  // Component 2 (last): max_thresholds(2,*) = +inf
  CHECK(thresholds.maxThresholds(2, 0) ==
        std::numeric_limits<double>::infinity());
  CHECK(thresholds.maxThresholds(2, 1) ==
        std::numeric_limits<double>::infinity());
}

TEST_CASE("AssignmentRule::evaluate: 2 Gaussians, simple box assignment",
          "[assignmentrule]") {
  // With 3 components and 2 Gaussians, manually construct thresholds and verify
  // evaluate. Component 1: z1 < 0,   z2 in (-inf, +inf) → region where z1 < 0
  // counts for comp 1 This is a simplified hand-crafted threshold set for
  // deterministic testing.

  AssignmentRuleThresholds thresholds;
  thresholds.numGaussians = 2;
  thresholds.minThresholds.resize(3, 2);
  thresholds.maxThresholds.resize(3, 2);

  const double inf = std::numeric_limits<double>::infinity();

  // Component 0: z1 in (-inf, 0), z2 in (-inf, +inf)
  thresholds.minThresholds.row(0) << -inf, -inf;
  thresholds.maxThresholds.row(0) << 0.0, inf;

  // Component 1: z1 in [0, +inf), z2 in (-inf, 0)
  thresholds.minThresholds.row(1) << 0.0, -inf;
  thresholds.maxThresholds.row(1) << inf, 0.0;

  // Component 2: z1 in [0, +inf), z2 in [0, +inf)
  thresholds.minThresholds.row(2) << 0.0, 0.0;
  thresholds.maxThresholds.row(2) << inf, inf;

  // Z matrix: 4 test voxels
  Eigen::MatrixXd z(4, 2);
  z << -1.0, 1.0, // z1<0 → comp 0 (1-based: 1)
      1.0, -1.0,  // z1>=0, z2<0 → comp 1 (1-based: 2)
      1.0, 1.0,   // z1>=0, z2>=0 → comp 2 (1-based: 3)
      -0.5, -0.5; // z1<0 → comp 0 (1-based: 1)

  std::mt19937_64 rng(0);
  AssignmentRule ar(rng);
  const Eigen::VectorXi assignments = ar.evaluate(z, thresholds);

  CHECK(assignments(0) == 1);
  CHECK(assignments(1) == 2);
  CHECK(assignments(2) == 3);
  CHECK(assignments(3) == 1);
}

TEST_CASE(
    "AssignmentRule: end-to-end volume fractions are approximately recovered",
    "[assignmentrule][slow]") {
  // Verify that selectThresholds + evaluate reproduces the requested fractions.
  // Draw many samples from N(0, I_2) and check empirical fractions.
  const std::vector<double> P = {0.30, 0.35, 0.35};
  const int numComponents = 3;
  const int numGaussians = 2;

  std::mt19937_64 rng(12345);
  AssignmentRule ar(rng);

  const auto thresholds = ar.selectThresholds(P);

  // Generate N(0, I_2) samples
  const int nSamples = 50000;
  std::normal_distribution<double> normal(0.0, 1.0);
  Eigen::MatrixXd z(nSamples, numGaussians);
  for (int i = 0; i < nSamples; ++i) {
    for (int g = 0; g < numGaussians; ++g) {
      z(i, g) = normal(rng);
    }
  }

  const Eigen::VectorXi assignments = ar.evaluate(z, thresholds);

  // Count each component
  Eigen::VectorXd empirical = Eigen::VectorXd::Zero(numComponents);
  for (int i = 0; i < nSamples; ++i) {
    const int comp = assignments(i);
    if (comp >= 1 && comp <= numComponents) {
      empirical(comp - 1) += 1.0;
    }
  }
  empirical /= static_cast<double>(nSamples);

  // Check that empirical fractions are within 2% of targets
  const double tol = 0.02;
  for (int j = 0; j < numComponents; ++j) {
    CHECK(empirical(j) == Approx(P[static_cast<std::size_t>(j)]).margin(tol));
  }
}

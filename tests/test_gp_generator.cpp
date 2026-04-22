#include "GPGenerator.hpp"
#include "PGRFSimulation.hpp"
#include "SimulationParams.hpp"

#include <Eigen/Dense>
#include <catch2/catch.hpp>
#include <cmath>
#include <numeric>
#include <random>

using namespace mtrsim;

// ─────────────────────────────────────────────────────────────────────────────
// GPGenerator tests
// ─────────────────────────────────────────────────────────────────────────────

namespace {
// Exponential correlation: rho(lag, theta) = exp(-|lag|/theta)
auto expCorrFn = [](double lag, double theta) -> double {
  return std::exp(-std::abs(lag) / theta);
};
} // anonymous namespace

TEST_CASE("GPGenerator: output length equals nx*ny*nz", "[gpgenerator]") {
  std::mt19937_64 rng(42);
  GPGenerator gpGen(rng, expCorrFn);

  auto gp = gpGen.generate(0.1, 0.1, 0.1, {1.0, 1.0, 1.0}, 10, 8, 3);
  CHECK(gp.size() == 10 * 8 * 3);
}

TEST_CASE("GPGenerator: 2D slab (nz=1) output length is nx*ny",
          "[gpgenerator]") {
  std::mt19937_64 rng(7);
  GPGenerator gpGen(rng, expCorrFn);

  auto gp = gpGen.generate(0.02, 0.02, 0.02, {0.1, 0.45, 0.1}, 20, 15, 1);
  CHECK(gp.size() == 20 * 15);
}

TEST_CASE("GPGenerator: same seed produces identical output", "[gpgenerator]") {
  auto gp1 = [&]() {
    std::mt19937_64 rng(999);
    GPGenerator gpGen(rng, expCorrFn);
    return gpGen.generate(0.1, 0.1, 0.1, {0.5, 0.5, 0.5}, 8, 8, 1);
  }();

  auto gp2 = [&]() {
    std::mt19937_64 rng(999);
    GPGenerator gpGen(rng, expCorrFn);
    return gpGen.generate(0.1, 0.1, 0.1, {0.5, 0.5, 0.5}, 8, 8, 1);
  }();

  CHECK(gp1.isApprox(gp2));
}

TEST_CASE("GPGenerator: different seeds produce different output",
          "[gpgenerator]") {
  std::mt19937_64 rng1(1), rng2(2);
  GPGenerator g1(rng1, expCorrFn), g2(rng2, expCorrFn);

  auto out1 = g1.generate(0.1, 0.1, 0.1, {0.5, 0.5, 0.5}, 10, 10, 1);
  auto out2 = g2.generate(0.1, 0.1, 0.1, {0.5, 0.5, 0.5}, 10, 10, 1);

  CHECK_FALSE(out1.isApprox(out2));
}

TEST_CASE("GPGenerator: sample mean ≈ 0 and variance ≈ 1 on large 1D grid",
          "[gpgenerator]") {
  // A large 1D grid (ny=nz=1) has enough voxels for reliable moment estimates.
  // The marginal distribution of each GP value is N(0,1) since Gamma(0,0) = 1.
  // Use a SHORT correlation length (theta=0.05) relative to spacing (h=0.02) so
  // the effective sample size is large: n_eff ≈ nx * (1-rho)/(1+rho) ≈ nx/6.
  std::mt19937_64 rng(12345);
  GPGenerator gpGen(rng, expCorrFn);

  const int nx = 5000;
  auto gp = gpGen.generate(0.02, 0.02, 0.02, {0.05, 0.05, 0.05}, nx, 1, 1);
  REQUIRE(gp.size() == nx);

  const double mean = gp.mean();
  const double var =
      (gp.array() - mean).square().sum() / static_cast<double>(nx - 1);

  // n_eff ≈ 5000/6 ≈ 833 → σ_mean ≈ 0.035, σ_var ≈ 0.05.  Use generous margins.
  CHECK(mean == Approx(0.0).margin(0.15));
  CHECK(var == Approx(1.0).margin(0.20));
}

TEST_CASE("GPGenerator: exponential autocorrelation approximately correct",
          "[gpgenerator]") {
  // For a 1D GP with exp covariance rho(h,theta) = exp(-h/theta),
  // the lag-1 autocorrelation should be ≈ exp(-spacing/theta).
  std::mt19937_64 rng(77777);
  GPGenerator gpGen(rng, expCorrFn);

  const int nx = 5000;
  const double spacing = 0.05;
  const double theta = 0.5;
  auto gp = gpGen.generate(spacing, spacing, spacing, {theta, theta, theta}, nx,
                           1, 1);

  // Estimate lag-1 autocorrelation
  const double mean = gp.mean();
  const Eigen::VectorXd centred = gp.array() - mean;
  double cov0 = centred.dot(centred) / static_cast<double>(nx);
  double cov1 = centred.head(nx - 1).dot(centred.tail(nx - 1)) /
                static_cast<double>(nx - 1);
  const double rho_estimated = cov1 / cov0;

  const double rho_expected = std::exp(-spacing / theta);
  CHECK(rho_estimated == Approx(rho_expected).margin(0.05));
}

// ─────────────────────────────────────────────────────────────────────────────
// PGRFSimulation tests
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("PGRFSimulation: output dimensions are correct", "[pgrfsimulation]") {
  SimulationParams params;
  // Small grid for speed
  params.xLen = 0.5;
  params.yLen = 0.5;
  params.zLen = 0.0;
  params.dx = 0.1;
  params.dy = 0.1;
  params.dz = 0.1;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.10, 0.10}, {0.10, 0.10, 0.10}};

  const int nx = static_cast<int>(std::round(params.xLen / params.dx)); // 5
  const int ny = static_cast<int>(std::round(params.yLen / params.dy)); // 5
  const int nz = 1;
  const int N = nx * ny * nz; // 25
  const int numGaussians =
      static_cast<int>(params.volumeFractions.size()) - 1; // 2

  std::mt19937_64 rng(42);
  PGRFSimulation sim(rng);
  const PGRFResult result = sim.run(params);

  CHECK(result.mtrIndex.size() == N);
  CHECK(result.latentFields.rows() == N);
  CHECK(result.latentFields.cols() == numGaussians);
}

TEST_CASE("PGRFSimulation: all voxel assignments are in valid range",
          "[pgrfsimulation]") {
  SimulationParams params;
  params.xLen = 0.5;
  params.yLen = 0.5;
  params.zLen = 0.0;
  params.dx = 0.1;
  params.dy = 0.1;
  params.dz = 0.1;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.10, 0.10}, {0.10, 0.10, 0.10}};

  const int numComponents = static_cast<int>(params.volumeFractions.size());

  std::mt19937_64 rng(99);
  PGRFSimulation sim(rng);
  const PGRFResult result = sim.run(params);

  CHECK((result.mtrIndex.array() >= 1).all());
  CHECK((result.mtrIndex.array() <= numComponents).all());
}

TEST_CASE("PGRFSimulation: volume fractions approximately match targets [slow]",
          "[pgrfsimulation][slow]") {
  // Use a grid large enough relative to the correlation length for reliable
  // estimates. theta=0.10, dx=0.02 → correlation range ≈ 3*theta = 0.3 mm. Grid
  // 4mm×4mm = 200×200 voxels → n_eff ≈ (4/0.3)^2 ≈ 178 independent patches.
  // σ(P̂1) ≈ sqrt(0.30*0.70/178) ≈ 0.034 → 10% margin is ~3 sigma.
  SimulationParams params;
  params.xLen = 4.0;
  params.yLen = 4.0;
  params.zLen = 0.0;
  params.dx = 0.02;
  params.dy = 0.02;
  params.dz = 0.02;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.10, 0.10}, {0.10, 0.10, 0.10}};

  const int nx = static_cast<int>(std::round(params.xLen / params.dx)); // 200
  const int ny = static_cast<int>(std::round(params.yLen / params.dy)); // 200
  const int N = nx * ny;
  const int numComponents = static_cast<int>(params.volumeFractions.size());

  std::mt19937_64 rng(42);
  PGRFSimulation sim(rng);
  const PGRFResult result = sim.run(params);

  // Allow generous tolerance: 3-sigma for the correlated field estimator
  const double tol = 0.10;
  for (int j = 1; j <= numComponents; ++j) {
    const int count = (result.mtrIndex.array() == j).count();
    const double empirical =
        static_cast<double>(count) / static_cast<double>(N);
    CHECK(empirical ==
          Approx(params.volumeFractions[static_cast<std::size_t>(j - 1)])
              .margin(tol));
  }
}

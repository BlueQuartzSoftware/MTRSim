#include "ISimulationObserver.hpp"
#include "MTRSimDriver.hpp"
#include "ODFCalculator.hpp"
#include "ODFSampler.hpp"

#include <Eigen/Dense>
#include <catch2/catch.hpp>
#include <cmath>
#include <numbers>
#include <numeric>
#include <random>

using namespace mtrsim;

// ─────────────────────────────────────────────────────────────────────────────
// CrystalSymmetry tests removed — orientation symmetry expansion now goes
// through EbsdLib's `LaueOps` / `Euler<double>` / `OrientationMatrix<double>`
// directly inside the consumers (ODFCalculator, ComputeODF). EbsdLib's own
// unit tests cover the expansion math.
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// ODFCalculator tests
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("ODFCalculator: output ODFComponent has correct sizes", "[odfcalculator]")
{
  ODFCalculator calc;

  // Small set of random orientations
  const int N = 10;
  Eigen::VectorXd phi1(N), phi(N), phi2(N);
  for(int i = 0; i < N; ++i)
  {
    phi1[i] = 0.1 * i;
    phi[i] = 0.05 * i;
    phi2[i] = 0.2 * i;
  }

  const ODFComponent odf = calc.compute(phi1, phi, phi2);

  const int expected = 72 * 36 * 72; // = 186624
  CHECK(odf.odfVal.size() == expected);
  CHECK(odf.phi1Bins.size() == expected);
  CHECK(odf.phiBins.size() == expected);
  CHECK(odf.phi2Bins.size() == expected);
}

TEST_CASE("ODFCalculator: odfVal is non-negative", "[odfcalculator]")
{
  ODFCalculator calc;

  Eigen::VectorXd phi1(3), phi(3), phi2(3);
  phi1 << 0.5, 1.0, 2.0;
  phi << 0.3, 0.6, 1.2;
  phi2 << 1.0, 2.0, 3.0;

  const ODFComponent odf = calc.compute(phi1, phi, phi2);

  CHECK((odf.odfVal.array() >= 0.0).all());
}

TEST_CASE("ODFCalculator: bin centres are within Euler-space range", "[odfcalculator]")
{
  ODFCalculator calc;

  Eigen::VectorXd phi1(2), phi(2), phi2(2);
  phi1 << 0.0, 1.0;
  phi << 0.0, 0.8;
  phi2 << 0.0, 2.0;

  const ODFComponent odf = calc.compute(phi1, phi, phi2);

  const double twoPi = 2.0 * std::numbers::pi;
  const double piVal = std::numbers::pi;

  CHECK((odf.phi1Bins.array() >= 0.0).all());
  CHECK((odf.phi1Bins.array() <= twoPi).all());
  CHECK((odf.phiBins.array() >= 0.0).all());
  CHECK((odf.phiBins.array() <= piVal).all());
  CHECK((odf.phi2Bins.array() >= 0.0).all());
  CHECK((odf.phi2Bins.array() <= twoPi).all());
}

TEST_CASE("ODFCalculator: sum of odfVal ≈ 1 for large grid", "[odfcalculator]")
{
  // The smoothing weights sum to 1.0 per orientation, so the total ODF mass
  // should equal 1.0 (matching MATLAB unnormalised output with count/N
  // scaling). Use a small number of orientations to keep the test fast.
  ODFCalculator calc;

  // 100 random-ish orientations distributed across Euler space
  const int N = 100;
  Eigen::VectorXd phi1(N), phi(N), phi2(N);
  for(int i = 0; i < N; ++i)
  {
    phi1[i] = 2.0 * std::numbers::pi * (static_cast<double>(i) / N);
    phi[i] = std::numbers::pi * (static_cast<double>(i) / N);
    phi2[i] = 2.0 * std::numbers::pi * (static_cast<double>((i * 37) % N) / N);
  }

  const ODFComponent odf = calc.compute(phi1, phi, phi2);

  // Total weight of smoothing kernel per orientation = 0.332 + 6*0.448/6 +
  // 12*0.16/12 + 8*0.06/8 = 1.0
  CHECK(odf.odfVal.sum() == Approx(1.0).margin(1e-9));
}

// ─────────────────────────────────────────────────────────────────────────────
// ODFSampler tests
// ─────────────────────────────────────────────────────────────────────────────

namespace
{
// Build a trivial ODFComponent: uniform probability across nBins bins,
// bin centres at (i+0.5)*binWidth in all three Euler directions.
ODFComponent makeUniformODF(int nBins, double binWidth)
{
  ODFComponent c;
  c.odfVal = Eigen::VectorXd::Ones(nBins) / static_cast<double>(nBins);
  c.phi1Bins = Eigen::VectorXd::LinSpaced(nBins, 0.5 * binWidth, (nBins - 0.5) * binWidth);
  c.phiBins = c.phi1Bins;
  c.phi2Bins = c.phi1Bins;
  return c;
}
} // anonymous namespace

TEST_CASE("ODFSampler::sampleN: output dimensions are N x 3", "[odfsampler]")
{
  std::mt19937_64 rng(42);
  ODFSampler sampler(rng);

  const int nBins = 100;
  const ODFComponent uODF = makeUniformODF(nBins, 0.05);

  const int N = 50;
  const Eigen::MatrixXd result = sampler.sampleN(N, uODF, uODF);

  CHECK(result.rows() == N);
  CHECK(result.cols() == 3);
}

TEST_CASE("ODFSampler::sampleN: same seed reproduces identical output", "[odfsampler]")
{
  const int nBins = 50;
  const ODFComponent uODF = makeUniformODF(nBins, 0.04);
  const int N = 20;

  Eigen::MatrixXd r1, r2;
  {
    std::mt19937_64 rng(7);
    ODFSampler sampler(rng);
    r1 = sampler.sampleN(N, uODF, uODF);
  }
  {
    std::mt19937_64 rng(7);
    ODFSampler sampler(rng);
    r2 = sampler.sampleN(N, uODF, uODF);
  }

  CHECK(r1.isApprox(r2));
}

TEST_CASE("ODFSampler::sampleN: uniform ODF samples cover bin centres "
          "uniformly [slow]",
          "[odfsampler][slow]")
{
  // With a uniform ODF and a large draw, each bin should appear roughly
  // equally.
  std::mt19937_64 rng(12345);
  ODFSampler sampler(rng);

  const int nBins = 10;
  const double bw = 2.0 * std::numbers::pi / nBins;
  const ODFComponent uODF = makeUniformODF(nBins, bw);

  const int N = 10000;
  const Eigen::MatrixXd result = sampler.sampleN(N, uODF, uODF);

  // Count how many samples fall in each bin by phi1
  Eigen::VectorXi counts = Eigen::VectorXi::Zero(nBins);
  for(int i = 0; i < N; ++i)
  {
    int bin = static_cast<int>(std::floor(result(i, 0) / bw));
    bin = std::clamp(bin, 0, nBins - 1);
    counts[bin]++;
  }

  // Each bin should receive ~N/nBins draws; allow ±30% relative tolerance
  const double expected = static_cast<double>(N) / nBins;
  for(int b = 0; b < nBins; ++b)
  {
    CHECK(static_cast<double>(counts[b]) == Approx(expected).margin(0.30 * expected));
  }
}

TEST_CASE("ODFSampler::sampleOne: returns a single valid orientation", "[odfsampler]")
{
  std::mt19937_64 rng(99);
  ODFSampler sampler(rng);

  const int nBins = 30;
  const double bw = std::numbers::pi / nBins;
  const ODFComponent uODF = makeUniformODF(nBins, bw);

  const EulerAngles ea = sampler.sampleOne(uODF, uODF);

  // Angles should be near the bin centres ± half a bin width
  CHECK(ea.phi1 >= -bw);
  CHECK(ea.phi >= -bw);
  CHECK(ea.phi2 >= -bw);
}

namespace
{
class ImmediateCancel : public mtrsim::ISimulationObserver
{
public:
  void updateProgress(int64_t, int64_t, const std::string&) override
  {
  }
  bool shouldCancel() const override
  {
    return true;
  }
};
} // namespace

TEST_CASE("sampleN bails out promptly when observer cancels", "[odf_sampler]")
{
  mtrsim::ODFComponent uni = mtrsim::buildUniformODF(72, 36, 72);
  std::mt19937_64 rng(1);
  mtrsim::ODFSampler sampler{rng};
  ImmediateCancel cancel;
  Eigen::MatrixXd out = sampler.sampleN(100000, uni, uni, &cancel);
  REQUIRE(out.rows() == 0);
}

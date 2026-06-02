#include "ISimulationObserver.hpp"
#include "MTRSimDriver.hpp"

#include <catch2/catch.hpp>

#include <array>
#include <numbers>

namespace
{
class CancelAfterObserver : public mtrsim::ISimulationObserver
{
public:
  explicit CancelAfterObserver(int k)
  : m_K(k)
  {
  }
  void updateProgress(int64_t, int64_t, const std::string&) override
  {
    ++m_Count;
  }
  bool shouldCancel() const override
  {
    return m_Count >= m_K;
  }

private:
  int m_K;
  mutable int m_Count = 0;
};
} // namespace

TEST_CASE("buildUniformODF produces correct bin centres", "[mtrsim_driver]")
{
  const mtrsim::ODFComponent uni = mtrsim::buildUniformODF(72, 36, 72);

  REQUIRE(uni.odfVal.size() == 72 * 36 * 72);
  REQUIRE(uni.phi1Bins.size() == 72 * 36 * 72);

  // Uniform mass: every bin equal, sums to 1.
  REQUIRE(uni.odfVal.sum() == Approx(1.0));
  REQUIRE(uni.odfVal[0] == Approx(1.0 / (72.0 * 36.0 * 72.0)));

  // First bin centre: i1=iPHI=i2=0 -> all 0.5 * step.
  REQUIRE(uni.phi1Bins[0] == Approx(0.5 * 2.0 * std::numbers::pi / 72.0));
  REQUIRE(uni.phiBins[0] == Approx(0.5 * std::numbers::pi / 36.0));
  REQUIRE(uni.phi2Bins[0] == Approx(0.5 * 2.0 * std::numbers::pi / 72.0));
}

TEST_CASE("gridToODFComponent derives bin centres in radians and normalizes", "[mtrsim_driver]")
{
  const int n1 = 72, nPHI = 36, n2 = 72;
  std::vector<double> values(n1 * nPHI * n2, 2.0); // unnormalized constant

  const mtrsim::ODFComponent c = mtrsim::gridToODFComponent(values, n1, nPHI, n2, 5.0, 5.0, 5.0);

  REQUIRE(c.odfVal.size() == n1 * nPHI * n2);
  REQUIRE(c.odfVal.sum() == Approx(1.0)); // normalized
  // 5 deg step -> first bin centre 2.5 deg in radians.
  const double deg2rad = std::numbers::pi / 180.0;
  REQUIRE(c.phi1Bins[0] == Approx(2.5 * deg2rad));
  REQUIRE(c.phiBins[0] == Approx(2.5 * deg2rad));
  REQUIRE(c.phi2Bins[0] == Approx(2.5 * deg2rad));
}

TEST_CASE("remapSimToZYX moves y-fastest data to x-fastest layout", "[mtrsim_driver]")
{
  // 2x3x2 grid (nx=2, ny=3, nz=2). Fill sim-order vector with its own index.
  const int nx = 2, ny = 3, nz = 2;
  std::vector<int> in(nx * ny * nz);
  for(int i = 0; i < nx * ny * nz; ++i)
  {
    in[i] = i;
  }

  const std::vector<int> out = mtrsim::remapSimToZYX(in, nx, ny, nz);

  // Spot-check: SIMPLNX (ix=1, iy=0, iz=0) -> kNx = 1.
  //   source sim index = iz*(nx*ny) + ix*ny + iy = 0 + 1*3 + 0 = 3.
  REQUIRE(out[1] == 3);
  // SIMPLNX (ix=0, iy=1, iz=0) -> kNx = 2; sim = 0 + 0 + 1 = 1.
  REQUIRE(out[2] == 1);
  // Non-zero iz slice: SIMPLNX (ix=1, iy=2, iz=1).
  //   kNx = 1*(3*2) + 2*2 + 1 = 6 + 4 + 1 = 11.
  //   sim = 1*(2*3) + 1*3 + 2 = 6 + 3 + 2 = 11.
  REQUIRE(out[11] == 11);
  // Non-zero iz slice: SIMPLNX (ix=0, iy=0, iz=1).
  //   kNx = 1*(3*2) + 0 + 0 = 6; sim = 1*(2*3) + 0 + 0 = 6.
  REQUIRE(out[6] == 6);
  // Same total size.
  REQUIRE(out.size() == in.size());
}

TEST_CASE("simulateMTR reproduces target volume fractions (statistical)", "[mtrsim_driver][statistical]")
{
  mtrsim::SimulationParams params;
  params.xLen = 6.0;
  params.yLen = 6.0;
  params.zLen = 0.0;
  params.dx = 0.02;
  params.dy = 0.02;
  params.dz = 0.02;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.45, 0.10}, {0.08, 0.37, 0.08}};
  params.seed = 42;

  std::vector<mtrsim::ODFComponent> comps = {mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72)};

  std::mt19937_64 rng(params.seed);
  const mtrsim::MTRSimResult r = mtrsim::simulateMTR(params, comps, rng, 72, 36, 72);

  const int N = r.nx * r.ny * r.nz;
  REQUIRE(static_cast<int>(r.mtrIndex.size()) == N);

  for(int v : r.mtrIndex)
  {
    REQUIRE(v >= 1);
    REQUIRE(v <= 3);
  }

  std::array<int, 3> counts{0, 0, 0};
  for(int v : r.mtrIndex)
  {
    counts[static_cast<std::size_t>(v - 1)]++;
  }
  REQUIRE(static_cast<double>(counts[0]) / N == Approx(0.30).margin(0.05));
  REQUIRE(static_cast<double>(counts[1]) / N == Approx(0.35).margin(0.05));
  REQUIRE(static_cast<double>(counts[2]) / N == Approx(0.35).margin(0.05));

  REQUIRE(static_cast<int>(r.phi1.size()) == N);
  REQUIRE(static_cast<int>(r.phi.size()) == N);
  REQUIRE(static_cast<int>(r.phi2.size()) == N);

  for(double a : r.phi1)
  {
    REQUIRE(a >= 0.0);
    REQUIRE(a <= 2.0 * std::numbers::pi);
  }
  for(double a : r.phi)
  {
    REQUIRE(a >= 0.0);
    REQUIRE(a <= std::numbers::pi);
  }
  for(double a : r.phi2)
  {
    REQUIRE(a >= 0.0);
    REQUIRE(a <= 2.0 * std::numbers::pi);
  }
}

TEST_CASE("simulateMTR cancels early when observer requests it", "[mtrsim_driver]")
{
  mtrsim::SimulationParams params;
  params.xLen = 6.0;
  params.yLen = 6.0;
  params.zLen = 0.0;
  params.dx = 0.02;
  params.dy = 0.02;
  params.dz = 0.02;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.45, 0.10}, {0.08, 0.37, 0.08}};
  params.seed = 42;
  std::vector<mtrsim::ODFComponent> comps = {mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72)};
  std::mt19937_64 rng(params.seed);
  CancelAfterObserver obs(1); // cancel at the first progress checkpoint
  const mtrsim::MTRSimResult r = mtrsim::simulateMTR(params, comps, rng, 72, 36, 72, &obs);
  REQUIRE(r.cancelled);
  REQUIRE(r.mtrIndex.empty());
}

TEST_CASE("simulateMTR with nullptr observer is unaffected", "[mtrsim_driver]")
{
  mtrsim::SimulationParams params;
  params.xLen = 2.0;
  params.yLen = 2.0;
  params.zLen = 0.0;
  params.dx = 0.02;
  params.dy = 0.02;
  params.dz = 0.02;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.45, 0.10}, {0.08, 0.37, 0.08}};
  std::vector<mtrsim::ODFComponent> comps = {mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72)};
  std::mt19937_64 rng(7);
  const mtrsim::MTRSimResult r = mtrsim::simulateMTR(params, comps, rng, 72, 36, 72);
  REQUIRE_FALSE(r.cancelled);
  REQUIRE_FALSE(r.mtrIndex.empty());
}

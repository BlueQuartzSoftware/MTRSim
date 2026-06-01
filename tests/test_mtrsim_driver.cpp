#include "MTRSimDriver.hpp"

#include <catch2/catch.hpp>

#include <numbers>

TEST_CASE("buildUniformODF produces correct bin centres", "[mtrsim_driver]") {
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

TEST_CASE("gridToODFComponent derives bin centres in radians and normalizes", "[mtrsim_driver]") {
  const int n1 = 72, nPHI = 36, n2 = 72;
  std::vector<double> values(n1 * nPHI * n2, 2.0); // unnormalized constant

  const mtrsim::ODFComponent c =
      mtrsim::gridToODFComponent(values, n1, nPHI, n2, 5.0, 5.0, 5.0);

  REQUIRE(c.odfVal.size() == n1 * nPHI * n2);
  REQUIRE(c.odfVal.sum() == Approx(1.0)); // normalized
  // 5 deg step -> first bin centre 2.5 deg in radians.
  const double deg2rad = std::numbers::pi / 180.0;
  REQUIRE(c.phi1Bins[0] == Approx(2.5 * deg2rad));
  REQUIRE(c.phiBins[0] == Approx(2.5 * deg2rad));
  REQUIRE(c.phi2Bins[0] == Approx(2.5 * deg2rad));
}

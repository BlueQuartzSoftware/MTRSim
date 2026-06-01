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

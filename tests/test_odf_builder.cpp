#include "ODFBuilder.hpp"

#include <catch2/catch.hpp>

#include <array>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

// -----------------------------------------------------------------------------
// Test 1: No-smoothing hits exactly one bin per tuple.
// -----------------------------------------------------------------------------

TEST_CASE(
    "ODFBuilder::accumulate (no smoothing) hits exactly one bin per tuple",
    "[ODFBuilder]") {
  constexpr double k_DegToRad = std::numbers::pi / 180.0;
  const mtrsim::ODFBuildParams p{72, 36, 72, 5.0, false};
  std::vector<double> values(72 * 36 * 72, 0.0);

  // Tuple at bin (i_phi1=5, i_PHI=3, i_phi2=10) center → angles
  // (27.5, 17.5, 52.5) deg
  const std::vector<std::array<double, 3>> eulers = {
      {27.5 * k_DegToRad, 17.5 * k_DegToRad, 52.5 * k_DegToRad}};
  mtrsim::accumulate(eulers, p, values);

  double total = 0.0;
  for (double v : values) {
    total += v;
  }
  REQUIRE(total == Approx(1.0).epsilon(1e-12));

  const std::size_t ix = static_cast<std::size_t>(5 * 36 * 72 + 3 * 72 + 10);
  REQUIRE(values[ix] == Approx(1.0).epsilon(1e-12));
}

// -----------------------------------------------------------------------------
// Test 2: Smoothing distributes per MATLAB weights, summing to 1.0.
// -----------------------------------------------------------------------------

TEST_CASE("ODFBuilder::accumulate (smoothing) distributes per MATLAB weights",
          "[ODFBuilder]") {
  constexpr double k_DegToRad = std::numbers::pi / 180.0;
  const mtrsim::ODFBuildParams p{72, 36, 72, 5.0, true};
  std::vector<double> values(72 * 36 * 72, 0.0);

  const std::vector<std::array<double, 3>> eulers = {
      {27.5 * k_DegToRad, 17.5 * k_DegToRad, 52.5 * k_DegToRad}};
  mtrsim::accumulate(eulers, p, values);

  // Total contribution per tuple = 0.332 + 0.448 + 0.16 + 0.06 = 1.0
  double total = 0.0;
  for (double v : values) {
    total += v;
  }
  REQUIRE(total == Approx(1.0).epsilon(1e-9));

  const std::size_t ix = static_cast<std::size_t>(5 * 36 * 72 + 3 * 72 + 10);
  REQUIRE(values[ix] == Approx(0.332).epsilon(1e-12));
}

// -----------------------------------------------------------------------------
// Test 3: normalize() divides in place.
// -----------------------------------------------------------------------------

TEST_CASE("ODFBuilder::normalize divides in place", "[ODFBuilder]") {
  std::vector<double> v = {2.0, 4.0, 8.0};
  mtrsim::normalize(v, 2.0);
  REQUIRE(v == std::vector<double>{1.0, 2.0, 4.0});
}

// -----------------------------------------------------------------------------
// Test 4: accumulate throws on values-size mismatch.
// -----------------------------------------------------------------------------

TEST_CASE("ODFBuilder::accumulate throws on values-size mismatch",
          "[ODFBuilder]") {
  const mtrsim::ODFBuildParams p{72, 36, 72, 5.0, false};
  std::vector<double> values(100, 0.0); // wrong size
  REQUIRE_THROWS_AS(mtrsim::accumulate({{0.1, 0.2, 0.3}}, p, values),
                    std::invalid_argument);
}

// -----------------------------------------------------------------------------
// Test 5: normalize with 0.0 is a no-op.
// -----------------------------------------------------------------------------

TEST_CASE("ODFBuilder::normalize is a no-op when normalizer is zero",
          "[ODFBuilder]") {
  std::vector<double> v = {1.0, 2.0, 3.0};
  mtrsim::normalize(v, 0.0);
  REQUIRE(v == std::vector<double>{1.0, 2.0, 3.0});
}

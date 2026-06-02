#include "ConfigIO.hpp"

#include <catch2/catch.hpp>

#include <cstdio>
#include <fstream>
#include <string>

namespace {
std::string writeTemp(const std::string& contents)
{
  const std::string path = std::string(MTRSIM_TEST_DATA_DIR) + "/_tmp_config_io.json";
  std::ofstream o(path);
  o << contents;
  o.close();
  return path;
}
} // namespace

TEST_CASE("parseConfigJson reads known fields and ignores extras", "[config_io]")
{
  const std::string path = writeTemp(R"({
    "xLen": 38.1, "yLen": 12.7, "zLen": 0.0,
    "dx": 0.02, "dy": 0.02, "dz": 0.02,
    "volumeFractions": [0.30, 0.35, 0.35],
    "thetaList": [[0.10,0.45,0.10],[0.08,0.37,0.08]],
    "odfInputPath": "ignored.h5", "nuggetVariance": [0.6,0.7,0.7], "seed": 99
  })");
  const mtrsim::SimulationParams p = mtrsim::parseConfigJson(path);
  REQUIRE(p.xLen == Approx(38.1));
  REQUIRE(p.dx == Approx(0.02));
  REQUIRE(p.volumeFractions.size() == 3);
  REQUIRE(p.volumeFractions[1] == Approx(0.35));
  REQUIRE(p.thetaList.size() == 2);
  REQUIRE(p.thetaList[0][1] == Approx(0.45));
  REQUIRE(p.seed == 99);
  std::remove(path.c_str());
}

TEST_CASE("parseConfigJson throws on missing file", "[config_io]")
{
  REQUIRE_THROWS_AS(mtrsim::parseConfigJson("/nonexistent/path/nope.json"), std::runtime_error);
}

TEST_CASE("parseConfigJson throws on malformed JSON", "[config_io]")
{
  const std::string path = writeTemp("{ this is not json");
  REQUIRE_THROWS_AS(mtrsim::parseConfigJson(path), std::runtime_error);
  std::remove(path.c_str());
}

#include "IPFMapper.hpp"
#include "ODFCalculator.hpp"
#include "PoleFigure.hpp"

#include <Eigen/Dense>
#include <catch2/catch.hpp>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <numbers>

using namespace mtrsim;

// ─────────────────────────────────────────────────────────────────────────────
// IPFMapper tests
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("IPFMapper::eulerToColors: output size equals input size", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  const int N = 15;
  Eigen::VectorXd phi1 = Eigen::VectorXd::LinSpaced(N, 0.0, 2.0 * std::numbers::pi);
  Eigen::VectorXd phi = Eigen::VectorXd::LinSpaced(N, 0.0, std::numbers::pi);
  Eigen::VectorXd phi2 = Eigen::VectorXd::LinSpaced(N, 0.0, 2.0 * std::numbers::pi);

  const auto colors = mapper.eulerToColors(phi1, phi, phi2);
  CHECK(static_cast<int>(colors.size()) == N);
}

TEST_CASE("IPFMapper::eulerToColors: all channel values in [0, 255]", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  const int N = 20;
  Eigen::VectorXd phi1(N), phi(N), phi2(N);
  for(int i = 0; i < N; ++i)
  {
    phi1[i] = 0.314159 * i;
    phi[i] = 0.157080 * i;
    phi2[i] = 0.628318 * i;
  }

  const auto colors = mapper.eulerToColors(phi1, phi, phi2);
  for(std::size_t k = 0; k < colors.size(); ++k)
  {
    // uint8_t is already [0,255] by type — this just forces the check to show up in test output
    CHECK(colors[k].r >= 0);
    CHECK(colors[k].g >= 0);
    CHECK(colors[k].b >= 0);
  }
}

TEST_CASE("IPFMapper::eulerToColors: identity orientation [0001] maps to red (EbsdLib)", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  Eigen::VectorXd phi1(1), phi(1), phi2(1);
  phi1[0] = 0.0;
  phi[0] = 0.0;
  phi2[0] = 1.0e-15;

  const auto colors = mapper.eulerToColors(phi1, phi, phi2);
  CHECK(colors[0].r == 255);
  CHECK(colors[0].g == 0);
  CHECK(colors[0].b == 0);
}

TEST_CASE("IPFMapper MatLab scheme: identity orientation [0001] maps to red", "[ipfmapper]")
{
  // phi1=0, PHI=0, phi2=0 → G=I → specimen normal [0,0,1] in crystal frame = c-axis.
  // Stereographic projection places [0001] at (X=0, Y=0), so r=0:
  //   cmap1 = 1, cmap2 = 0, cmap3 = 0  →  RGB = (255, 0, 0).
  IPFMapper mapper{CrystalSystem::HCP};

  Eigen::VectorXd phi1(1), phi(1), phi2(1);
  phi1[0] = 0.0;
  phi[0] = 0.0;
  phi2[0] = 1.0e-15;

  const auto colors = mapper.eulerToColors(phi1, phi, phi2, {0.0, 0.0, 1.0}, IPFColorScheme::MatLab);
  CHECK(colors[0].r == 255);
  CHECK(colors[0].g == 0);
  CHECK(colors[0].b == 0);
}

TEST_CASE("IPFMapper MatLab scheme: output size equals input size", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  const int N = 15;
  Eigen::VectorXd phi1 = Eigen::VectorXd::LinSpaced(N, 0.0, 2.0 * std::numbers::pi);
  Eigen::VectorXd phi = Eigen::VectorXd::LinSpaced(N, 0.0, std::numbers::pi);
  Eigen::VectorXd phi2 = Eigen::VectorXd::LinSpaced(N, 0.0, 2.0 * std::numbers::pi);

  const auto colors = mapper.eulerToColors(phi1, phi, phi2, {0.0, 0.0, 1.0}, IPFColorScheme::MatLab);
  CHECK(static_cast<int>(colors.size()) == N);
}

TEST_CASE("IPFMapper MatLab scheme: same orientation gives identical colours", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  const double p1 = 0.7;
  const double ph = 0.4;
  const double p2 = 1.1;

  Eigen::VectorXd phi1(3), phi(3), phi2(3);
  phi1 << p1, p1, p1;
  phi << ph, ph, ph;
  phi2 << p2, p2, p2;

  const auto colors = mapper.eulerToColors(phi1, phi, phi2, {0.0, 0.0, 1.0}, IPFColorScheme::MatLab);
  CHECK(colors[0].r == colors[1].r);
  CHECK(colors[0].g == colors[1].g);
  CHECK(colors[0].b == colors[1].b);
  CHECK(colors[0].r == colors[2].r);
  CHECK(colors[0].g == colors[2].g);
  CHECK(colors[0].b == colors[2].b);
}

TEST_CASE("IPFMapper MatLab scheme: FCC throws", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::FCC};

  Eigen::VectorXd phi1(1), phi(1), phi2(1);
  phi1[0] = 0.0;
  phi[0] = 0.0;
  phi2[0] = 0.1;

  CHECK_THROWS_AS(mapper.eulerToColors(phi1, phi, phi2, {0.0, 0.0, 1.0}, IPFColorScheme::MatLab), std::invalid_argument);
}

TEST_CASE("IPFMapper::writeIPFTriangleLegendMatLab: creates file on disk", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};
  const std::string outPath = "/tmp/test_ipf_triangle_legend.png";

  REQUIRE_NOTHROW(mapper.writeIPFTriangleLegendMatLab(256, outPath));
  CHECK(std::filesystem::exists(outPath));

  // Verify it has non-zero size
  CHECK(std::filesystem::file_size(outPath) > 0);
  std::filesystem::remove(outPath);
}

TEST_CASE("IPFMapper::writeIPFTriangleLegendMatLab: FCC throws", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::FCC};
  CHECK_THROWS_AS(mapper.writeIPFTriangleLegendMatLab(256, "/tmp/should_not_exist.png"), std::invalid_argument);
  CHECK_FALSE(std::filesystem::exists("/tmp/should_not_exist.png"));
}

TEST_CASE("IPFMapper: compare EbsdLib vs MatLab colours", "[ipfmapper][.print]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  // Representative Euler angles (radians) spanning different orientations
  struct EulerSet
  {
    double phi1, phi, phi2;
    const char* label;
  };

  // clang-format off
  const std::array<EulerSet, 8> angles = {{
    {0.0,   0.0,   1e-15, "[0001] c-axis (identity)"},
    {0.0,   1.571, 0.0,   "PHI=90deg"},
    {0.5,   0.3,   0.8,   "arbitrary #1"},
    {1.2,   0.7,   2.1,   "arbitrary #2"},
    {2.0,   1.0,   0.5,   "arbitrary #3"},
    {0.0,   0.5,   1.047, "phi2=60deg"},
    {3.14,  1.57,  0.0,   "near 180/90/0"},
    {0.785, 0.524, 0.524, "45/30/30 deg"},
  }};
  // clang-format on

  std::printf("\n%-30s | %-17s | %-17s\n", "Orientation", "EbsdLib (R,G,B)", "MatLab  (R,G,B)");
  std::printf("-------------------------------+-------------------+------------------\n");

  for(const auto& a : angles)
  {
    Eigen::VectorXd p1(1), p(1), p2(1);
    p1[0] = a.phi1;
    p[0] = a.phi;
    p2[0] = a.phi2;

    const auto ebsd = mapper.eulerToColors(p1, p, p2, {0.0, 0.0, 1.0}, IPFColorScheme::EbsdLib);
    const auto matlab = mapper.eulerToColors(p1, p, p2, {0.0, 0.0, 1.0}, IPFColorScheme::MatLab);

    std::printf("%-30s | (%3d, %3d, %3d)   | (%3d, %3d, %3d)\n", a.label, ebsd[0].r, ebsd[0].g, ebsd[0].b, matlab[0].r, matlab[0].g, matlab[0].b);
  }

  std::printf("\n");
  CHECK(true); // keep Catch2 happy
}

TEST_CASE("IPFMapper: [0001] to [2-1-10] sweep (PHI 0-90 by 1 deg)", "[ipfmapper][.print]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  // Walking from [0001] to [2-1-10]:
  //   phi1 = 0, phi2 = 0, PHI steps from 0° to 90°.
  // At PHI=0 the specimen normal aligns with [0001] → red.
  // At PHI=90° it aligns with [2-1-10] → green.
  constexpr double k_Deg2Rad = std::numbers::pi / 180.0;

  std::printf("\n PHI (deg) | EbsdLib (R,G,B)   | MatLab  (R,G,B)\n");
  std::printf("-----------+-------------------+------------------\n");

  for(int deg = 0; deg <= 90; ++deg)
  {
    Eigen::VectorXd p1(1), p(1), p2(1);
    p1[0] = 0.0;
    p[0] = deg * k_Deg2Rad;
    p2[0] = (deg == 0) ? 1.0e-15 : 0.0; // nudge at exactly 0 to match existing convention

    const auto ebsd = mapper.eulerToColors(p1, p, p2, {0.0, 0.0, 1.0}, IPFColorScheme::EbsdLib);
    const auto matlab = mapper.eulerToColors(p1, p, p2, {0.0, 0.0, 1.0}, IPFColorScheme::MatLab);

    std::printf("    %3d    | (%3d, %3d, %3d)   | (%3d, %3d, %3d)\n", deg, ebsd[0].r, ebsd[0].g, ebsd[0].b, matlab[0].r, matlab[0].g, matlab[0].b);
  }

  std::printf("\n");
  CHECK(true);
}

TEST_CASE("IPFMapper::eulerToColors: same orientation gives identical colours", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  const double p1 = 0.7;
  const double ph = 0.4;
  const double p2 = 1.1;

  Eigen::VectorXd phi1(3), phi(3), phi2(3);
  phi1 << p1, p1, p1;
  phi << ph, ph, ph;
  phi2 << p2, p2, p2;

  const auto colors = mapper.eulerToColors(phi1, phi, phi2);
  CHECK(colors[0].r == colors[1].r);
  CHECK(colors[0].g == colors[1].g);
  CHECK(colors[0].b == colors[1].b);
  CHECK(colors[0].r == colors[2].r);
  CHECK(colors[0].g == colors[2].g);
  CHECK(colors[0].b == colors[2].b);
}

TEST_CASE("IPFMapper::eulerToColors: size mismatch throws", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  Eigen::VectorXd phi1(3), phi(2), phi2(3);
  phi1.setZero();
  phi.setZero();
  phi2.setZero();

  CHECK_THROWS_AS(mapper.eulerToColors(phi1, phi, phi2), std::invalid_argument);
}

TEST_CASE("IPFMapper::writePNG: creates file on disk", "[ipfmapper]")
{
  IPFMapper mapper{CrystalSystem::HCP};

  // 3×3 regular spatial grid
  const int nx = 3, ny = 3;
  const int N = nx * ny;

  Eigen::MatrixXd coords(N, 2);
  Eigen::VectorXd phi1(N), phi(N), phi2(N);
  int idx = 0;
  for(int iy = 0; iy < ny; ++iy)
  {
    for(int ix = 0; ix < nx; ++ix)
    {
      coords(idx, 0) = static_cast<double>(ix) * 0.1;
      coords(idx, 1) = static_cast<double>(iy) * 0.1;
      phi1[idx] = 0.2 * idx;
      phi[idx] = 0.1 * idx;
      phi2[idx] = 0.3 * idx;
      ++idx;
    }
  }

  const std::string outPath = "/tmp/test_ipf_mapper_output.png";
  REQUIRE_NOTHROW(mapper.writePNG(coords, phi1, phi, phi2, outPath));
  CHECK(std::filesystem::exists(outPath));
  std::filesystem::remove(outPath);
}

// ─────────────────────────────────────────────────────────────────────────────
// PoleFigure tests
// ─────────────────────────────────────────────────────────────────────────────

namespace
{
// Build a small ODFComponent with a few non-zero bins.
ODFComponent makeTinyODF()
{
  ODFCalculator calc;
  Eigen::VectorXd phi1(5), phi(5), phi2(5);
  phi1 << 0.3, 1.0, 2.0, 3.5, 5.0;
  phi << 0.2, 0.5, 0.8, 1.1, 1.4;
  phi2 << 0.6, 1.5, 2.5, 3.0, 4.5;
  return calc.compute(phi1, phi, phi2);
}
} // anonymous namespace

TEST_CASE("PoleFigure::fromODF: returns non-empty data", "[polefigure]")
{
  PoleFigure pf;
  const ODFComponent odf = makeTinyODF();
  const PoleFigureData pfd = pf.fromODF(odf);

  CHECK(pfd.x.size() > 0);
  CHECK(pfd.y.size() == pfd.x.size());
  CHECK(pfd.intensity.size() == pfd.x.size());
}

TEST_CASE("PoleFigure::fromODF: all intensities non-negative", "[polefigure]")
{
  PoleFigure pf;
  const ODFComponent odf = makeTinyODF();
  const PoleFigureData pfd = pf.fromODF(odf);

  CHECK((pfd.intensity.array() >= 0.0).all());
}

TEST_CASE("PoleFigure::fromODF: intensities normalised to sum ~= 1", "[polefigure]")
{
  PoleFigure pf;
  const ODFComponent odf = makeTinyODF();
  const PoleFigureData pfd = pf.fromODF(odf);

  CHECK(pfd.intensity.sum() == Approx(1.0).margin(1.0e-9));
}

TEST_CASE("PoleFigure::fromODF: stereographic coords are finite", "[polefigure]")
{
  PoleFigure pf;
  const ODFComponent odf = makeTinyODF();
  const PoleFigureData pfd = pf.fromODF(odf);

  // All projected X, Y values should be finite (no NaN/Inf)
  for(int i = 0; i < pfd.x.size(); ++i)
  {
    CHECK(std::isfinite(pfd.x[i]));
    CHECK(std::isfinite(pfd.y[i]));
  }
}

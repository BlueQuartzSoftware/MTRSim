// Port of unit_triangle_IPF_coords.m, IPF_colors.m, and view_IPF_map.m.
// Daniel M. Sparkman, 07/05/2017.

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "IPFMapper.hpp"

#include <EbsdLib/LaueOps/CubicOps.h>
#include <EbsdLib/LaueOps/HexagonalOps.h>
#include <EbsdLib/Utilities/ColorTable.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>

#if MTRSIM_HAS_TBB
#include <tbb/parallel_for.h>
#endif

namespace mtrsim
{

// ─────────────────────────────────────────────────────────────────────────────
// Symmetry-operator Euler-angle tables (radians) — from MATLAB source.
// Each row is { phi1, PHI, phi2 }.
// ─────────────────────────────────────────────────────────────────────────────

namespace
{
constexpr double k_Deg2Rad = std::numbers::pi / 180.0;

// HCP — 12 operators
constexpr std::array<std::array<double, 3>, 12> k_HcpSymOps = {{
    {0, 0, 0},
    {0, 0, 60.0 * k_Deg2Rad},
    {0, 0, 120.0 * k_Deg2Rad},
    {0, 0, 180.0 * k_Deg2Rad},
    {0, 0, 240.0 * k_Deg2Rad},
    {0, 0, 300.0 * k_Deg2Rad},
    {0, 180.0 * k_Deg2Rad, 0},
    {0, 180.0 * k_Deg2Rad, 60.0 * k_Deg2Rad},
    {0, 180.0 * k_Deg2Rad, 120.0 * k_Deg2Rad},
    {0, 180.0 * k_Deg2Rad, 180.0 * k_Deg2Rad},
    {0, 180.0 * k_Deg2Rad, 240.0 * k_Deg2Rad},
    {0, 180.0 * k_Deg2Rad, 300.0 * k_Deg2Rad},
}};

// ─────────────────────────────────────────────────────────────────────────────
// Build a Bunge passive-rotation matrix from Euler angles (phi1, PHI, phi2).
//
//   G = | c1c2−s1s2C    s1c2+c1s2C    s2S  |
//       | −c1s2−s1c2C  −s1s2+c1c2C   c2S  |
//       |  s1S          −c1S          C    |
// ─────────────────────────────────────────────────────────────────────────────
inline Eigen::Matrix3d bungeRotationMatrix(double p1, double P, double p2)
{
  const double c1 = std::cos(p1), s1 = std::sin(p1);
  const double C = std::cos(P), S = std::sin(P);
  const double c2 = std::cos(p2), s2 = std::sin(p2);

  Eigen::Matrix3d G;
  G(0, 0) = c1 * c2 - s1 * s2 * C;
  G(0, 1) = s1 * c2 + c1 * s2 * C;
  G(0, 2) = s2 * S;
  G(1, 0) = -c1 * s2 - s1 * c2 * C;
  G(1, 1) = -s1 * s2 + c1 * c2 * C;
  G(1, 2) = c2 * S;
  G(2, 0) = s1 * S;
  G(2, 1) = -c1 * S;
  G(2, 2) = C;
  return G;
}
} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
IPFMapper::IPFMapper(CrystalSystem system)
: m_System(system)
{
}

// ─────────────────────────────────────────────────────────────────────────────
// Public dispatch — routes to the selected colour scheme.

std::vector<RGBColor> IPFMapper::eulerToColors(const Eigen::VectorXd& phi1, const Eigen::VectorXd& phi, const Eigen::VectorXd& phi2, std::array<double, 3> refDir, IPFColorScheme scheme) const
{
  const int N = static_cast<int>(phi1.size());
  if(phi.size() != N || phi2.size() != N)
  {
    throw std::invalid_argument("IPFMapper::eulerToColors — phi1, phi, phi2 must have the same length");
  }

  if(scheme == IPFColorScheme::MatLab)
  {
    return eulerToColorsMatLab(phi1, phi, phi2);
  }
  return eulerToColorsEbsdLib(phi1, phi, phi2, refDir);
}

// ─────────────────────────────────────────────────────────────────────────────
// EbsdLib colour path — delegates to LaueOps::generateIPFColor().

std::vector<RGBColor> IPFMapper::eulerToColorsEbsdLib(const Eigen::VectorXd& phi1, const Eigen::VectorXd& phi, const Eigen::VectorXd& phi2, std::array<double, 3> refDir) const
{
  const int N = static_cast<int>(phi1.size());

  std::shared_ptr<ebsdlib::LaueOps> ops;
  if(m_System == CrystalSystem::HCP)
  {
    ops = ebsdlib::HexagonalOps::New();
  }
  else if(m_System == CrystalSystem::FCC)
  {
    ops = ebsdlib::CubicOps::New();
  }
  else
  {
    throw std::invalid_argument("IPFMapper::eulerToColorsEbsdLib — unsupported crystal system");
  }

  std::vector<RGBColor> colors(static_cast<std::size_t>(N));

  auto computeColor = [&](int i) {
    double eulers[3] = {phi1[i], phi[i], phi2[i]};
    const ebsdlib::Rgb argb = ops->generateIPFColor(eulers, refDir.data(), false);
    const std::size_t idx = static_cast<std::size_t>(i);
    colors[idx].r = static_cast<uint8_t>(ebsdlib::RgbColor::dRed(argb));
    colors[idx].g = static_cast<uint8_t>(ebsdlib::RgbColor::dGreen(argb));
    colors[idx].b = static_cast<uint8_t>(ebsdlib::RgbColor::dBlue(argb));
  };

#if MTRSIM_HAS_TBB
  tbb::parallel_for(0, N, computeColor);
#else
  for(int i = 0; i < N; ++i)
  {
    computeColor(i);
  }
#endif

  return colors;
}

// ─────────────────────────────────────────────────────────────────────────────
// MATLAB colour path — port of unit_triangle_IPF_coords.m + IPF_colors.m.
//
// Algorithm:
//   1. Build the Bunge rotation matrix G from each orientation.
//   2. For each crystal symmetry operator O_k, compute g_rot = O_k * G.
//   3. h_rot = g_rot * [0,0,1]  (third column of g_rot).
//   4. Flip to lower hemisphere if h_rot.z > 0.
//   5. Stereographic projection:  X = hx/(1−hz),  Y = hy/(1−hz).
//   6. Pick the first operator whose (X,Y) lands in the fundamental zone:
//        X ≥ 0,  Y ≥ 0,  Y ≤ X·tan(π/6)    [HCP]
//   7. Polar colour mapping (IPF_colors.m):
//        r     = √(X² + Y²)
//        β     = atan2(Y, X)
//        red   = (1 − r)  /  max(…)
//        green = r·(1 − β/(π/6))  /  max(…)
//        blue  = r·β/(π/6)  /  max(…)

std::vector<RGBColor> IPFMapper::eulerToColorsMatLab(const Eigen::VectorXd& phi1, const Eigen::VectorXd& phi, const Eigen::VectorXd& phi2) const
{
  if(m_System != CrystalSystem::HCP)
  {
    throw std::invalid_argument("IPFMapper::eulerToColorsMatLab — MatLab colour scheme currently supports HCP only");
  }

  const int N = static_cast<int>(phi1.size());
  constexpr int numOps = static_cast<int>(k_HcpSymOps.size());
  constexpr double k_MaxAngle = std::numbers::pi / 6.0; // 30° — HCP fundamental-zone arc
  const double tanMaxAngle = std::tan(k_MaxAngle);

  // Pre-compute symmetry-operator rotation matrices
  std::array<Eigen::Matrix3d, 12> symMats;
  for(int k = 0; k < numOps; ++k)
  {
    symMats[static_cast<std::size_t>(k)] = bungeRotationMatrix(k_HcpSymOps[static_cast<std::size_t>(k)][0], k_HcpSymOps[static_cast<std::size_t>(k)][1], k_HcpSymOps[static_cast<std::size_t>(k)][2]);
  }

  std::vector<RGBColor> colors(static_cast<std::size_t>(N));

  auto computeColor = [&](int i) {
    // Replicate the MATLAB guard: nudge phi2 away from exact zero
    double p2 = phi2[i];
    if(p2 == 0.0)
    {
      p2 += 1.0e-15;
    }

    const Eigen::Matrix3d G = bungeRotationMatrix(phi1[i], phi[i], p2);

    double xFund = 0.0;
    double yFund = 0.0;

    for(int k = 0; k < numOps; ++k)
    {
      // g_rot_symm = O_k * G   (matches MATLAB exactly)
      const Eigen::Matrix3d R = symMats[static_cast<std::size_t>(k)] * G;

      // h_rot = R * [0,0,1] = third column of R
      double hx = R(0, 2);
      double hy = R(1, 2);
      double hz = R(2, 2);

      // Flip to lower hemisphere (MATLAB: h_rot_3 > 0 → negate)
      if(hz > 0.0)
      {
        hx = -hx;
        hy = -hy;
        hz = -hz;
      }

      // Stereographic projection
      const double denom = 1.0 - hz; // always ≥ 1.0 after hemisphere flip
      const double X = hx / denom;
      const double Y = hy / denom;

      // Fundamental-zone check (HCP unit triangle)
      if(X >= 0.0 && Y >= 0.0 && Y <= X * tanMaxAngle)
      {
        xFund = X;
        yFund = Y;
        break;
      }
    }

    // ── Polar colour mapping (IPF_colors.m) ──────────────────────────────────
    const double r = std::sqrt(xFund * xFund + yFund * yFund);
    const double beta = std::atan2(yFund, xFund);

    double cmap1 = 1.0 - r;
    double cmap2 = r * (1.0 - beta / k_MaxAngle);
    double cmap3 = r * (beta / k_MaxAngle);

    // Normalise so that the brightest channel is 1.0
    const double maxC = std::max({cmap1, cmap2, cmap3});
    if(maxC > 0.0)
    {
      cmap1 /= maxC;
      cmap2 /= maxC;
      cmap3 /= maxC;
    }

    const std::size_t idx = static_cast<std::size_t>(i);
    colors[idx].r = static_cast<uint8_t>(std::clamp(std::lround(cmap1 * 255.0), 0L, 255L));
    colors[idx].g = static_cast<uint8_t>(std::clamp(std::lround(cmap2 * 255.0), 0L, 255L));
    colors[idx].b = static_cast<uint8_t>(std::clamp(std::lround(cmap3 * 255.0), 0L, 255L));
  };

#if MTRSIM_HAS_TBB
  tbb::parallel_for(0, N, computeColor);
#else
  for(int i = 0; i < N; ++i)
  {
    computeColor(i);
  }
#endif

  return colors;
}

// ─────────────────────────────────────────────────────────────────────────────
// writePNG — renders the IPF colour map to a PNG file.

void IPFMapper::writePNG(const Eigen::MatrixXd& spatialCoords, const Eigen::VectorXd& phi1_in, const Eigen::VectorXd& phi_in, const Eigen::VectorXd& phi2_in, const std::string& outputPath,
                         IPFColorScheme scheme) const
{
  const int N = static_cast<int>(phi1_in.size());
  if(spatialCoords.rows() != N || spatialCoords.cols() < 2)
  {
    throw std::invalid_argument("IPFMapper::writePNG — spatialCoords must have N rows and ≥ 2 columns");
  }

  const std::vector<RGBColor> colors = eulerToColors(phi1_in, phi_in, phi2_in, {0.0, 0.0, 1.0}, scheme);

  // ── Build sorted unique coordinate lists ────────────────────────────────────
  const Eigen::VectorXd xCol = spatialCoords.col(0);
  const Eigen::VectorXd yCol = spatialCoords.col(1);

  const double xRange = xCol.maxCoeff() - xCol.minCoeff();
  const double yRange = yCol.maxCoeff() - yCol.minCoeff();
  const double tol = 1.0e-6 * std::max({xRange, yRange, 1.0});

  auto sortedUnique = [&](const Eigen::VectorXd& v) -> std::vector<double> {
    std::vector<double> vals(v.data(), v.data() + v.size());
    std::sort(vals.begin(), vals.end());
    std::vector<double> unique;
    for(const double val : vals)
    {
      if(unique.empty() || std::abs(val - unique.back()) > tol)
      {
        unique.push_back(val);
      }
    }
    return unique;
  };

  const std::vector<double> xUnique = sortedUnique(xCol);
  const std::vector<double> yUnique = sortedUnique(yCol);
  const int width = static_cast<int>(xUnique.size());
  const int height = static_cast<int>(yUnique.size());

  if(width < 1 || height < 1)
  {
    throw std::runtime_error("IPFMapper::writePNG — could not determine image dimensions from spatialCoords");
  }

  // ── Fill pixel buffer (RGB, 3 bytes per pixel) ──────────────────────────────
  std::vector<uint8_t> pixels(static_cast<std::size_t>(width * height * 3), 0u);

  auto findRank = [&](const std::vector<double>& sorted, double val) -> int {
    const auto it = std::lower_bound(sorted.begin(), sorted.end(), val - tol);
    return static_cast<int>(std::distance(sorted.begin(), it));
  };

  for(int i = 0; i < N; ++i)
  {
    const int col = findRank(xUnique, xCol[i]);
    const int row = findRank(yUnique, yCol[i]);
    if(col < 0 || col >= width || row < 0 || row >= height)
    {
      continue;
    }
    const std::size_t px = static_cast<std::size_t>((row * width + col) * 3);
    const std::size_t ci = static_cast<std::size_t>(i);
    pixels[px] = colors[ci].r;
    pixels[px + 1] = colors[ci].g;
    pixels[px + 2] = colors[ci].b;
  }

  // ── Write PNG ───────────────────────────────────────────────────────────────
  const int stride = width * 3;
  if(stbi_write_png(outputPath.c_str(), width, height, 3, pixels.data(), stride) == 0)
  {
    throw std::runtime_error("IPFMapper::writePNG — stbi_write_png failed: " + outputPath);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// writeIPFTriangleLegendMatLab — renders the HCP fundamental-zone triangle
// using the MATLAB polar colour mapping (IPF_colors.m) and writes it as PNG.
//
// The HCP stereographic unit triangle is bounded by:
//   • X ≥ 0             (left edge:  [0001] → [10-10])
//   • Y ≥ 0             (bottom edge: [0001] → [2-1-10])
//   • Y ≤ X·tan(π/6)   (upper edge:  [0001] → [10-10])
//   • X² + Y² ≤ 1      (arc:         [2-1-10] → [10-10])

void IPFMapper::writeIPFTriangleLegendMatLab(int imageDim, const std::string& outputPath) const
{
  if(m_System != CrystalSystem::HCP)
  {
    throw std::invalid_argument("IPFMapper::writeIPFTriangleLegendMatLab — currently supports HCP only");
  }

  constexpr double k_MaxAngle = std::numbers::pi / 6.0;
  const double tanMaxAngle = std::tan(k_MaxAngle);
  const double sinMaxAngle = std::sin(k_MaxAngle);

  // Coordinate bounds with a small margin so the triangle edges aren't clipped.
  constexpr double k_Pad = 0.05;
  constexpr double xMin = -k_Pad;
  constexpr double xMax = 1.0 + k_Pad;
  const double yMin = -k_Pad;
  const double yMax = sinMaxAngle + k_Pad;
  const double xRange = xMax - xMin;
  const double yRange = yMax - yMin;

  const int width = imageDim;
  const int height = std::max(1, static_cast<int>(std::round(imageDim * yRange / xRange)));

  // White background
  std::vector<uint8_t> pixels(static_cast<std::size_t>(width * height * 3), 255u);

  for(int row = 0; row < height; ++row)
  {
    for(int col = 0; col < width; ++col)
    {
      // Map pixel centre to stereographic coordinates (Y increases upward)
      const double X = xMin + (static_cast<double>(col) + 0.5) * xRange / width;
      const double Y = yMax - (static_cast<double>(row) + 0.5) * yRange / height;

      // Fundamental-zone test
      if(X < 0.0 || Y < 0.0 || Y > X * tanMaxAngle || (X * X + Y * Y) > 1.0)
      {
        continue;
      }

      // Polar colour mapping — identical to IPF_colors.m
      const double r = std::sqrt(X * X + Y * Y);
      const double beta = std::atan2(Y, X);

      double cmap1 = 1.0 - r;
      double cmap2 = r * (1.0 - beta / k_MaxAngle);
      double cmap3 = r * (beta / k_MaxAngle);

      const double maxC = std::max({cmap1, cmap2, cmap3});
      if(maxC > 0.0)
      {
        cmap1 /= maxC;
        cmap2 /= maxC;
        cmap3 /= maxC;
      }

      const std::size_t px = static_cast<std::size_t>((row * width + col) * 3);
      pixels[px] = static_cast<uint8_t>(std::clamp(std::lround(cmap1 * 255.0), 0L, 255L));
      pixels[px + 1] = static_cast<uint8_t>(std::clamp(std::lround(cmap2 * 255.0), 0L, 255L));
      pixels[px + 2] = static_cast<uint8_t>(std::clamp(std::lround(cmap3 * 255.0), 0L, 255L));
    }
  }

  const int stride = width * 3;
  if(stbi_write_png(outputPath.c_str(), width, height, 3, pixels.data(), stride) == 0)
  {
    throw std::runtime_error("IPFMapper::writeIPFTriangleLegendMatLab — stbi_write_png failed: " + outputPath);
  }
}

} // namespace mtrsim

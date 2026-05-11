#include "ODFBuilder.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>

namespace mtrsim
{

namespace
{

// -----------------------------------------------------------------------------
// MATLAB calc_ODF.m tri-linear smoothing weights.
// Sum: 0.332 + 0.448 + 0.16 + 0.06 = 1.0 (exact in IEEE-754 double).
// -----------------------------------------------------------------------------
constexpr double k_CenterWeight = 0.332;
constexpr double k_FaceWeight = 0.448 / 6.0;
constexpr double k_EdgeWeight = 0.16 / 12.0;
constexpr double k_CornerWeight = 0.06 / 8.0;

// 6 face neighbors: +/-1 along exactly one axis.
constexpr int k_FaceOffsets[6][3] = {{-1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, +1, 0}, {0, 0, -1}, {0, 0, +1}};

// 12 edge neighbors: +/-1 along exactly two axes.
constexpr int k_EdgeOffsets[12][3] = {
    {-1, -1, 0}, {-1, +1, 0}, {+1, -1, 0}, {+1, +1, 0}, {-1, 0, -1}, {-1, 0, +1}, {+1, 0, -1}, {+1, 0, +1}, {0, -1, -1}, {0, -1, +1}, {0, +1, -1}, {0, +1, +1}};

/// Modulo wrap that handles negative dividends. Used unchanged for the
/// +1/-1 smoothing-neighbor stencil on all three Bunge axes -- this matches
/// calc_ODF.m's jf_minus/jf_plus/kf_minus/kf_plus/lf_minus/lf_plus logic
/// (lines 95-113), which IS uniform modulo on every axis.
inline int32_t wrap(int32_t i, int32_t n)
{
  return ((i % n) + n) % n;
}

/// Clamp a bin index to [0, n-1]. Matches MATLAB calc_ODF.m's upper-bound
/// clamping for bin assignment (lines 43-48): a sample exactly at the upper
/// axis bound (2*pi for phi1/phi2, pi for PHI) goes to the LAST bin, not
/// wrap-around to bin 0. Negative inputs (which shouldn't occur for
/// non-negative Bunge angles) defensively wrap via wrap().
inline int32_t clampBin(int32_t i, int32_t n)
{
  if(i < 0)
  {
    return wrap(i, n);
  }
  if(i >= n)
  {
    return n - 1;
  }
  return i;
}

/// Row-major linearization: i_phi1 * (nPHI * nphi2) + i_PHI * nphi2 + i_phi2.
inline std::size_t linearize(int32_t iPhi1, int32_t iPHI, int32_t iPhi2, int32_t nPHI, int32_t nphi2)
{
  return static_cast<std::size_t>(iPhi1) * static_cast<std::size_t>(nPHI) * static_cast<std::size_t>(nphi2) + static_cast<std::size_t>(iPHI) * static_cast<std::size_t>(nphi2)
         + static_cast<std::size_t>(iPhi2);
}

} // namespace

void accumulate(const std::vector<std::array<double, 3>>& eulersRad, const ODFBuildParams& params, std::vector<double>& values)
{
  const std::size_t expectedSize = static_cast<std::size_t>(params.nphi1) * static_cast<std::size_t>(params.nPHI) * static_cast<std::size_t>(params.nphi2);
  if(values.size() != expectedSize)
  {
    throw std::invalid_argument("ODFBuilder::accumulate: values size mismatch");
  }

  constexpr double k_RadToDeg = 180.0 / std::numbers::pi;

  for(const auto& tuple : eulersRad)
  {
    const double phi1Deg = tuple[0] * k_RadToDeg;
    const double PHIDeg = tuple[1] * k_RadToDeg;
    const double phi2Deg = tuple[2] * k_RadToDeg;

    // MATLAB calc_ODF.m clamps at the upper bound for bin assignment (lines 43-48):
    // an angle exactly equal to 2*pi (phi1/phi2) or pi (PHI) goes to the LAST bin,
    // not wrap-around to bin 0. Smoothing-neighbor identification (below) still
    // uses uniform modulo wrap on all three axes, matching calc_ODF.m's
    // jf_minus/jf_plus/kf_minus/kf_plus/lf_minus/lf_plus logic (lines 95-113).
    const int32_t iPhi1 = clampBin(static_cast<int32_t>(std::floor(phi1Deg / params.binSizeDeg)), params.nphi1);
    const int32_t iPHI = clampBin(static_cast<int32_t>(std::floor(PHIDeg / params.binSizeDeg)), params.nPHI);
    const int32_t iPhi2 = clampBin(static_cast<int32_t>(std::floor(phi2Deg / params.binSizeDeg)), params.nphi2);

    if(!params.smoothing)
    {
      values[linearize(iPhi1, iPHI, iPhi2, params.nPHI, params.nphi2)] += 1.0;
      continue;
    }

    // Center.
    values[linearize(iPhi1, iPHI, iPhi2, params.nPHI, params.nphi2)] += k_CenterWeight;

    // Faces (6).
    for(const auto& offs : k_FaceOffsets)
    {
      const int32_t jj = wrap(iPhi1 + offs[0], params.nphi1);
      const int32_t kk = wrap(iPHI + offs[1], params.nPHI);
      const int32_t ll = wrap(iPhi2 + offs[2], params.nphi2);
      values[linearize(jj, kk, ll, params.nPHI, params.nphi2)] += k_FaceWeight;
    }

    // Edges (12).
    for(const auto& offs : k_EdgeOffsets)
    {
      const int32_t jj = wrap(iPhi1 + offs[0], params.nphi1);
      const int32_t kk = wrap(iPHI + offs[1], params.nPHI);
      const int32_t ll = wrap(iPhi2 + offs[2], params.nphi2);
      values[linearize(jj, kk, ll, params.nPHI, params.nphi2)] += k_EdgeWeight;
    }

    // Corners (8): +/-1 along all three axes.
    for(int dj = -1; dj <= 1; dj += 2)
    {
      for(int dk = -1; dk <= 1; dk += 2)
      {
        for(int dl = -1; dl <= 1; dl += 2)
        {
          const int32_t jj = wrap(iPhi1 + dj, params.nphi1);
          const int32_t kk = wrap(iPHI + dk, params.nPHI);
          const int32_t ll = wrap(iPhi2 + dl, params.nphi2);
          values[linearize(jj, kk, ll, params.nPHI, params.nphi2)] += k_CornerWeight;
        }
      }
    }
  }
}

void normalize(std::vector<double>& values, double normalizer)
{
  if(normalizer == 0.0)
  {
    return;
  }
  for(auto& v : values)
  {
    v /= normalizer;
  }
}

} // namespace mtrsim

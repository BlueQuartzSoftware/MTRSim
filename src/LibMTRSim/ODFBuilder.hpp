#pragma once

#include "libmtrsim_export.h"

#include <array>
#include <cstdint>
#include <vector>

namespace mtrsim
{

/**
 * @brief Discretization parameters used by ODFBuilder::accumulate.
 *
 * Describes a uniform Bunge-Euler histogram grid with row-major layout
 *   linear = i_phi1 * (nPHI * nphi2) + i_PHI * nphi2 + i_phi2
 * (phi1 is the slowest-varying axis, phi2 the fastest).
 */
struct LIBMTRSIM_EXPORT ODFBuildParams
{
  int32_t nphi1;     ///< Number of bins along phi1 (slowest-varying, Z in ImageGeom).
  int32_t nPHI;      ///< Number of bins along PHI  (middle,           Y).
  int32_t nphi2;     ///< Number of bins along phi2 (fastest-varying,  X).
  double binSizeDeg; ///< Uniform bin size in degrees (all three axes).
  bool smoothing;    ///< When true, distribute each tuple over 27 bins (tri-linear
                     ///< smoothing).
};

/**
 * @brief Accumulates each Euler tuple's unit contribution into @p values.
 *
 * Inputs are Bunge Euler angles in RADIANS. The @p values buffer must be
 * pre-sized to nphi1 * nPHI * nphi2 (row-major, phi1-slowest) and initially
 * zero — this function only accumulates, it never clears. It does NOT
 * normalize.
 *
 * When smoothing is disabled, each tuple contributes 1.0 to its center bin.
 * When smoothing is enabled, the unit contribution is distributed across 27
 * bins (1 center + 6 faces + 12 edges + 8 corners) using the MATLAB
 * calc_ODF.m weights: 0.332 / 0.448/6 / 0.16/12 / 0.06/8 (summing to 1.0).
 *
 * @throws std::invalid_argument if @p values.size() does not match
 *         nphi1 * nPHI * nphi2.
 */
LIBMTRSIM_EXPORT void accumulate(const std::vector<std::array<double, 3>>& eulersRad, const ODFBuildParams& params, std::vector<double>& values);

/**
 * @brief In-place division: values[i] /= normalizer.
 *
 * No-op when @p normalizer equals 0.0 (avoids divide-by-zero blow-ups from
 * empty accumulators).
 */
LIBMTRSIM_EXPORT void normalize(std::vector<double>& values, double normalizer);

} // namespace mtrsim

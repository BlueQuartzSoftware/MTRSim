// Port of sample_N_orientations_from_ODF.m and sample_orientation_from_ODF.m.
// Daniel M. Sparkman, Research.

#include "ODFSampler.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace mtrsim {

ODFSampler::ODFSampler(std::mt19937_64 &rng) : m_Rng(rng) {}

// ─────────────────────────────────────────────────────────────────────────────
// sampleN — port of sample_N_orientations_from_ODF.m.
//
// Algorithm:
//   1. Build CDF from component.odfVal (normalised to sum = 1).
//   2. Draw N uniform [0,1] samples; find bin via upper_bound on the CDF.
//   3. Shuffle the bin-index array (matching MATLAB randperm).
//   4. For each draw, take the bin centre from uniform.phi1Bins / phiBins /
//      phi2Bins and add independent jitter ∈ [-binWidth/2, +binWidth/2].
//
// The bin width is inferred as uniform.phi1Bins[1] − uniform.phi1Bins[0],
// matching MATLAB's bin_scaling_factor = degree_spacing/180*pi.

Eigen::MatrixXd ODFSampler::sampleN(int n, const ODFComponent &component,
                                    const ODFComponent &uniform) {
  const int nBins = static_cast<int>(component.odfVal.size());
  if (nBins < 1) {
    throw std::invalid_argument("ODFSampler::sampleN — odfVal is empty");
  }
  if (uniform.phi1Bins.size() != static_cast<Eigen::Index>(nBins) ||
      uniform.phiBins.size() != static_cast<Eigen::Index>(nBins) ||
      uniform.phi2Bins.size() != static_cast<Eigen::Index>(nBins)) {
    throw std::invalid_argument(
        "ODFSampler::sampleN — uniform bin vectors do not match odfVal size");
  }

  // ── Build normalised CDF ──────────────────────────────────────────────────
  const double total = component.odfVal.sum();
  std::vector<double> cdf(static_cast<std::size_t>(nBins) + 1);
  cdf[0] = 0.0;
  for (int j = 0; j < nBins; ++j) {
    cdf[static_cast<std::size_t>(j) + 1] =
        cdf[static_cast<std::size_t>(j)] + component.odfVal[j] / total;
  }

  // ── Draw N uniform samples and map to bin indices ─────────────────────────
  std::uniform_real_distribution<double> uDist(0.0, 1.0);
  std::vector<int> binIdx(static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    const double u = uDist(m_Rng);
    // upper_bound finds the first position where cdf > u;
    // the bin index is one before that position.
    const auto it = std::upper_bound(cdf.cbegin(), cdf.cend(), u);
    int bin = static_cast<int>(std::distance(cdf.cbegin(), it)) - 1;
    binIdx[static_cast<std::size_t>(i)] = std::clamp(bin, 0, nBins - 1);
  }

  // ── Shuffle (matches MATLAB randperm) ─────────────────────────────────────
  std::shuffle(binIdx.begin(), binIdx.end(), m_Rng);

  // ── Infer bin width from the uniform ODF bin centres ─────────────────────
  // For standard 5° spacing: binWidth ≈ 5*π/180 rad.
  const double binWidth = (nBins > 1)
                              ? (uniform.phi1Bins[1] - uniform.phi1Bins[0])
                              : (5.0 * std::acos(-1.0) / 180.0);

  // ── Assign bin-centre coordinates with independent per-angle jitter ───────
  std::uniform_real_distribution<double> jitter(-0.5 * binWidth,
                                                0.5 * binWidth);

  Eigen::MatrixXd result(n, 3);
  for (int i = 0; i < n; ++i) {
    const int b = binIdx[static_cast<std::size_t>(i)];
    result(i, 0) = uniform.phi1Bins[b] + jitter(m_Rng);
    result(i, 1) = uniform.phiBins[b] + jitter(m_Rng);
    result(i, 2) = uniform.phi2Bins[b] + jitter(m_Rng);
  }

  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// sampleOne — port of sample_orientation_from_ODF.m.
//
// Uses the same inverse-CDF approach as sampleN but draws a single orientation.

EulerAngles ODFSampler::sampleOne(const ODFComponent &component,
                                  const ODFComponent &uniform) {
  const Eigen::MatrixXd row = sampleN(1, component, uniform);
  return EulerAngles{row(0, 0), row(0, 1), row(0, 2)};
}

} // namespace mtrsim

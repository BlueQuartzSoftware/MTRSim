// Port of calc_ODF.m — smoothed 3-D Euler-space histogram ODF calculation.
// Daniel M. Sparkman, 08/13/2013.
//
// All orientation math goes through EbsdLib's `Euler<double>`,
// `OrientationMatrix<double>`, and `LaueOps` directly. No local wrappers.
// All math is double-precision regardless of input precision (matching
// EbsdLib's internal convention).

#include "ODFCalculator.hpp"

#include <EbsdLib/Core/EbsdLibConstants.h>
#include <EbsdLib/LaueOps/LaueOps.h>
#include <EbsdLib/Math/Matrix3X3.hpp>
#include <EbsdLib/Orientation/Euler.hpp>
#include <EbsdLib/Orientation/OrientationFwd.hpp>
#include <EbsdLib/Orientation/OrientationMatrix.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace mtrsim
{

// ─────────────────────────────────────────────────────────────────────────────
// compute — port of calc_ODF.m (always-executed vectorised branch,
// smooth_ODF=true).
//
// Grid layout (matching MATLAB, degree_spacing = 5 by default):
//   phi1  : 72 bins covering [0, 2π)
//   PHI   : 36 bins covering [0, π)
//   phi2  : 72 bins covering [0, 2π)
//   Total : 72 × 36 × 72 = 186 624 bins
//
// Flat index (0-based, same ordering as MATLAB 1-based formula minus 1):
//   ix = phi1_ix * (nPHI * nphi2) + PHI_ix * nphi2 + phi2_ix
//
// Smoothing weights: centre 0.332, face 0.448/6, edge 0.16/12, corner 0.06/8.
// Wrapping: periodic in phi1 and phi2; PHI also wraps for the smoothing
// neighbours (matching MATLAB behaviour).

ODFComponent ODFCalculator::compute(const Eigen::VectorXd& phi1in, const Eigen::VectorXd& phiIn, const Eigen::VectorXd& phi2in, double degSpacing)
{
  if(phi1in.size() != phiIn.size() || phi1in.size() != phi2in.size())
  {
    throw std::invalid_argument("ODFCalculator::compute — phi1, phi, phi2 must have the same length");
  }

  // ── Grid constants ────────────────────────────────────────────────────────
  const int nBins1 = static_cast<int>(std::round(360.0 / degSpacing)); // 72
  const int nBinsPHI = nBins1 / 2;                                     // 36
  const int nBins2 = nBins1;                                           // 72
  const int nTotal = nBins1 * nBinsPHI * nBins2;                       // 186 624

  const double radSpacing = degSpacing * std::numbers::pi / 180.0;
  const double k_TwoPiOver = 2.0 * std::numbers::pi / static_cast<double>(nBins1);
  const double k_PiOver = std::numbers::pi / static_cast<double>(nBinsPHI);

  // ── Pre-compute bin centres for the output ODFComponent ──────────────────
  Eigen::VectorXd phi1Bins(nTotal);
  Eigen::VectorXd phiBins(nTotal);
  Eigen::VectorXd phi2Bins(nTotal);

  for(int ix = 0; ix < nTotal; ++ix)
  {
    const int i1 = ix / (nBinsPHI * nBins2);
    const int iPHI = (ix % (nBinsPHI * nBins2)) / nBins2;
    const int i2 = ix % nBins2;
    phi1Bins[ix] = (i1 + 0.5) * k_TwoPiOver;
    phiBins[ix] = (iPHI + 0.5) * k_PiOver;
    phi2Bins[ix] = (i2 + 0.5) * k_TwoPiOver;
  }

  // ── Symmetry expansion via EbsdLib (Hexagonal_High = MATLAB default) ──────
  // Idiom mirrors `simplnx`'s OrientationAnalysis filters (e.g. ComputeGBCD,
  // ComputeGBCDMetricBased): construct an `EulerD`, transform to an
  // `OrientationMatrix<double>`, and compose with `LaueOps::getMatSymOpD(k)`.
  const std::vector<ebsdlib::LaueOps::Pointer> allOps = ebsdlib::LaueOps::GetAllOrientationOps();
  const ebsdlib::LaueOps::Pointer laueOps = allOps[ebsdlib::CrystalStructure::Hexagonal_High];
  const std::size_t numOps = laueOps->getNumSymOps();

  const int N_input = static_cast<int>(phi1in.size());
  const int N = N_input * static_cast<int>(numOps);

  // ── Smoothing factors ─────────────────────────────────────────────────────
  const double normFactor = 1.0 / static_cast<double>(N);
  const double centerFactor = 0.332;
  const double faceFactor = 0.448 / 6.0;
  const double edgeFactor = 0.16 / 12.0;
  const double cornerFactor = 0.06 / 8.0;

  // ── Helper lambdas for periodic neighbour wrapping ────────────────────────
  // phi1 and phi2 are periodic over nBins1/nBins2 bins.
  // PHI also wraps for the smoothing kernel (matching MATLAB).
  auto wrap1 = [nBins1](int idx) -> int { return (idx % nBins1 + nBins1) % nBins1; };
  auto wrapPHI = [nBinsPHI](int idx) -> int { return (idx % nBinsPHI + nBinsPHI) % nBinsPHI; };
  auto wrap2 = [nBins2](int idx) -> int { return (idx % nBins2 + nBins2) % nBins2; };

  auto flatIdx = [nBinsPHI, nBins2](int i1, int iPHI, int i2) -> int { return i1 * nBinsPHI * nBins2 + iPHI * nBins2 + i2; };

  // ── Accumulate ODF ────────────────────────────────────────────────────────
  Eigen::VectorXd odfVal = Eigen::VectorXd::Zero(nTotal);

  for(int i = 0; i < N_input; ++i)
  {
    // Build the passive G matrix for the input orientation. EbsdLib's
    // `Euler::toOrientationMatrix()` applies the |x|<1e-7 → 0 eps-snap
    // we want to match the MATLAB reference.
    const ebsdlib::EulerDType eu(phi1in[i], phiIn[i], phi2in[i]);
    const ebsdlib::Matrix3X3D Gpassive = eu.toOrientationMatrix().toGMatrix();

    for(std::size_t k = 0; k < numOps; ++k)
    {
      // `getMatSymOpD(k)` returns the ACTIVE form; transpose for passive.
      // Symmetric variant: R = O_passive * G_passive_input.
      const ebsdlib::Matrix3X3D Op_active = laueOps->getMatSymOpD(k);
      const ebsdlib::Matrix3X3D R = Op_active.transpose() * Gpassive;

      // Extract Euler back via canonical om2eu (handles degenerate-PHI
      // and the eps-snap correctly, matching MATLAB after commit 8bbd6e).
      const ebsdlib::OrientationMatrixDType om(R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
      const ebsdlib::EulerDType outEu = om.toEuler();

      // Compute 0-based bin indices (clamp to valid range)
      const int jf = std::min(static_cast<int>(std::trunc(outEu[0] / radSpacing)), nBins1 - 1);
      const int kf = std::min(static_cast<int>(std::trunc(outEu[1] / radSpacing)), nBinsPHI - 1);
      const int lf = std::min(static_cast<int>(std::trunc(outEu[2] / radSpacing)), nBins2 - 1);

      // Wrapped neighbour indices
      const int jm = wrap1(jf - 1);
      const int jp = wrap1(jf + 1);
      const int km = wrapPHI(kf - 1);
      const int kp = wrapPHI(kf + 1);
      const int lm = wrap2(lf - 1);
      const int lp = wrap2(lf + 1);

      const double w = normFactor;

      // Centre bin
      odfVal(flatIdx(jf, kf, lf)) += w * centerFactor;

      // 6 face bins
      odfVal(flatIdx(jm, kf, lf)) += w * faceFactor;
      odfVal(flatIdx(jp, kf, lf)) += w * faceFactor;
      odfVal(flatIdx(jf, km, lf)) += w * faceFactor;
      odfVal(flatIdx(jf, kp, lf)) += w * faceFactor;
      odfVal(flatIdx(jf, kf, lm)) += w * faceFactor;
      odfVal(flatIdx(jf, kf, lp)) += w * faceFactor;

      // 12 edge bins
      odfVal(flatIdx(jm, kf, lm)) += w * edgeFactor;
      odfVal(flatIdx(jm, km, lf)) += w * edgeFactor;
      odfVal(flatIdx(jm, kp, lf)) += w * edgeFactor;
      odfVal(flatIdx(jm, kf, lp)) += w * edgeFactor;
      odfVal(flatIdx(jf, km, lm)) += w * edgeFactor;
      odfVal(flatIdx(jf, km, lp)) += w * edgeFactor;
      odfVal(flatIdx(jf, kp, lp)) += w * edgeFactor;
      odfVal(flatIdx(jf, kp, lm)) += w * edgeFactor;
      odfVal(flatIdx(jp, kf, lm)) += w * edgeFactor;
      odfVal(flatIdx(jp, km, lf)) += w * edgeFactor;
      odfVal(flatIdx(jp, kp, lf)) += w * edgeFactor;
      odfVal(flatIdx(jp, kf, lp)) += w * edgeFactor;

      // 8 corner bins
      odfVal(flatIdx(jm, km, lm)) += w * cornerFactor;
      odfVal(flatIdx(jm, km, lp)) += w * cornerFactor;
      odfVal(flatIdx(jm, kp, lp)) += w * cornerFactor;
      odfVal(flatIdx(jm, kp, lm)) += w * cornerFactor;
      odfVal(flatIdx(jp, km, lm)) += w * cornerFactor;
      odfVal(flatIdx(jp, km, lp)) += w * cornerFactor;
      odfVal(flatIdx(jp, kp, lp)) += w * cornerFactor;
      odfVal(flatIdx(jp, kp, lm)) += w * cornerFactor;
    }
  }

  ODFComponent result;
  result.odfVal = std::move(odfVal);
  result.phi1Bins = std::move(phi1Bins);
  result.phiBins = std::move(phiBins);
  result.phi2Bins = std::move(phi2Bins);
  return result;
}

} // namespace mtrsim

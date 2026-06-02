// Port of convert_ODF_to_PF.m — compute [0001] pole figure from an ODF.
// Daniel M. Sparkman, 08/18/2013.

#include "PoleFigure.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace mtrsim
{

// ─────────────────────────────────────────────────────────────────────────────
// fromODF — port of convert_ODF_to_PF.m (direct computation, no prebuilt
// tables).
//
// For each non-zero ODF bin (phi1, PHI, phi2) with value v:
//   1. Compute the [0001] c-axis direction in sample frame:
//        h_sample = G^T * [0,0,1] = last row of G = [sin(phi1)*sin(PHI),
//        −cos(phi1)*sin(PHI), cos(PHI)]
//   2. Stereographic projection from south pole (centre = [0001] up, PHI = 0):
//        X = h_x / (1 + h_z),  Y = h_y / (1 + h_z)
//   3. Intensity weight (Euler-space → pole-figure Jacobian, matching MATLAB
//   formula):
//        weight = v * 4π² / (dphi1 * dphi2 * |cos(PHI − dPHI/2) − cos(PHI +
//        dPHI/2)|)
//   4. Collect all (X, Y, weight) pairs; normalise intensities to sum = 1.
//
// The degSpacing parameter must match the spacing used to compute the ODF.

PoleFigureData PoleFigure::fromODF(const ODFComponent& component, double degSpacing)
{
  const int nTotal = static_cast<int>(component.odfVal.size());
  if(nTotal < 1)
  {
    throw std::invalid_argument("PoleFigure::fromODF — empty ODFComponent");
  }
  if(component.phi1Bins.size() != static_cast<Eigen::Index>(nTotal) || component.phiBins.size() != static_cast<Eigen::Index>(nTotal) || component.phi2Bins.size() != static_cast<Eigen::Index>(nTotal))
  {
    throw std::invalid_argument("PoleFigure::fromODF — ODFComponent bin "
                                "vectors have inconsistent sizes");
  }

  // ── Euler-space bin spacings
  // ───────────────────────────────────────────────── dphi1 = dphi2 = dPHI =
  // degSpacing * π / 180 for the standard layout.
  const double radSpacing = degSpacing * std::numbers::pi / 180.0;
  const double dphi1 = radSpacing;
  const double dphi2 = radSpacing;
  const double dPHI = radSpacing;
  const double halfDPHI = 0.5 * dPHI;

  // Prefactor from MATLAB: 4π² / (dphi1 * dphi2)
  const double prefactor = 4.0 * std::numbers::pi * std::numbers::pi / (dphi1 * dphi2);

  // ── Accumulate (X, Y, weight) for all non-zero bins
  // ──────────────────────────
  std::vector<double> xVec;
  std::vector<double> yVec;
  std::vector<double> wVec;
  xVec.reserve(static_cast<std::size_t>(nTotal / 4));
  yVec.reserve(xVec.capacity());
  wVec.reserve(xVec.capacity());

  double totalWeight = 0.0;

  for(int m = 0; m < nTotal; ++m)
  {
    const double v = component.odfVal[m];
    if(v <= 0.0)
    {
      continue;
    }

    const double phi1_m = component.phi1Bins[m];
    const double PHI_m = component.phiBins[m];

    // ── [0001] direction in sample frame: last row of G ──────────────────────
    // G[2,0] = sin(phi1)*sin(PHI), G[2,1] = −cos(phi1)*sin(PHI), G[2,2] =
    // cos(PHI)
    const double sinPHI = std::sin(PHI_m);
    const double cosPHI = std::cos(PHI_m);
    const double h_x = std::sin(phi1_m) * sinPHI;
    const double h_y = -std::cos(phi1_m) * sinPHI;
    const double h_z = cosPHI;

    // ── Stereographic projection from south pole
    // ────────────────────────────── Singularity when h_z → −1 (PHI → π): skip
    // those bins.
    const double denom = 1.0 + h_z;
    if(std::abs(denom) < 1.0e-8)
    {
      continue;
    }
    const double X = h_x / denom;
    const double Y = h_y / denom;

    // ── Jacobian weight (MATLAB convert_ODF_to_PF.m formula) ─────────────────
    // |cos(PHI − dPHI/2) − cos(PHI + dPHI/2)| = 2·|sin(PHI)|·sin(dPHI/2)
    const double cosDiff = std::abs(std::cos(PHI_m - halfDPHI) - std::cos(PHI_m + halfDPHI));
    if(cosDiff < 1.0e-12)
    {
      continue; // PHI ≈ 0 or π: vanishing angular area, skip
    }
    const double weight = v * prefactor / cosDiff;

    xVec.push_back(X);
    yVec.push_back(Y);
    wVec.push_back(weight);
    totalWeight += weight;
  }

  // ── Normalise intensities
  // ─────────────────────────────────────────────────────
  if(totalWeight > 0.0)
  {
    for(double& w : wVec)
    {
      w /= totalWeight;
    }
  }

  // ── Pack into PoleFigureData
  // ──────────────────────────────────────────────────
  const int nOut = static_cast<int>(xVec.size());
  PoleFigureData pf;
  pf.x = Eigen::Map<Eigen::VectorXd>(xVec.data(), nOut);
  pf.y = Eigen::Map<Eigen::VectorXd>(yVec.data(), nOut);
  pf.intensity = Eigen::Map<Eigen::VectorXd>(wVec.data(), nOut);
  return pf;
}

} // namespace mtrsim

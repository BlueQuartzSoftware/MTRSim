// Port of symmetric_euler_angles.m — crystallographic symmetry expansion.
// Daniel M. Sparkman, 02/14/2013.

#include "CrystalSymmetry.hpp"

#include <EbsdLib/LaueOps/CubicOps.h>
#include <EbsdLib/LaueOps/HexagonalOps.h>
#include <EbsdLib/Math/Matrix3X3.hpp>

#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>

namespace mtrsim
{

CrystalSymmetry::CrystalSymmetry(CrystalSystem system)
: m_System(system)
{
}

int CrystalSymmetry::numOperators() const
{
  switch(m_System)
  {
  case CrystalSystem::HCP:
    return 12;
  case CrystalSystem::FCC:
    return 24;
  default:
    return 0;
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// expand — port of symmetric_euler_angles.m (all_symmetry = true branch).
//
// For each input orientation (phi1, PHI, phi2) and each symmetry operator k,
// computes the composite passive rotation matrix R = O_k * G, then extracts
// Euler angles from R.
//
// EbsdLib convention: LaueOps::getMatSymOpD(k) returns the ACTIVE (alibi)
// rotation matrix, which is the transpose of the passive rotation matrix used
// in the MATLAB Bunge convention.  So O_passive = getMatSymOpD(k).transpose().
//
// The Bunge eu2om formula (matching EbsdLib's OrientationTransformation::eu2om):
//
//   G = | c1c2-s1Cs2    s1c2+c1Cs2    Ss2 |
//       | -c1s2-s1Cc2  -s1s2+c1Cc2   Sc2 |
//       |  s1S          -c1S           C  |
//
// Euler angle extraction from R (EbsdLib's OrientationTransformation::om2eu):
//   PHI  = acos(R[2,2])
//   phi1 = atan2(R[2,0], -R[2,1])
//   phi2 = atan2(R[0,2],  R[1,2])

Eigen::MatrixXd CrystalSymmetry::expand(const Eigen::VectorXd& phi1, const Eigen::VectorXd& phi, const Eigen::VectorXd& phi2) const
{
  const int N = static_cast<int>(phi1.size());
  const int numOps = numOperators();

  // Instantiate the appropriate LaueOps object
  std::shared_ptr<ebsdlib::LaueOps> ops;
  if(m_System == CrystalSystem::HCP)
  {
    ops = ebsdlib::HexagonalOps::New();
  }
  else
  {
    ops = ebsdlib::CubicOps::New();
  }

  Eigen::MatrixXd result(N * numOps, 3);

  constexpr double k_Eps = 1.0e-6;
  const double k_Pi = std::numbers::pi;
  const double k_TwoPi = 2.0 * k_Pi;

  for(int i = 0; i < N; ++i)
  {
    // ── Build passive rotation matrix G from Bunge Euler angles ─────────────
    const double c1 = std::cos(phi1[i]);
    const double s1 = std::sin(phi1[i]);
    const double C = std::cos(phi[i]);
    const double S = std::sin(phi[i]);
    const double c2 = std::cos(phi2[i]);
    const double s2 = std::sin(phi2[i]);

    // Row-major layout: G[row*3 + col]
    const ebsdlib::Matrix3X3D G(c1 * c2 - s1 * s2 * C, s1 * c2 + c1 * s2 * C, s2 * S, -c1 * s2 - s1 * c2 * C, -s1 * s2 + c1 * c2 * C, c2 * S, s1 * S, -c1 * S, C);

    for(int k = 0; k < numOps; ++k)
    {
      // getMatSymOpD(k) stores the ACTIVE (= O_passive^T) rotation matrix.
      // O_passive = O_t.transpose()
      // R = O_passive * G = O_t.transpose() * G
      const ebsdlib::Matrix3X3D O_t = ops->getMatSymOpD(static_cast<size_t>(k));
      const ebsdlib::Matrix3X3D R = O_t.transpose() * G;

      // ── Extract Euler angles from R (EbsdLib's om2eu convention) ───────────
      // Indices in row-major storage: R[8]=R(2,2), R[6]=R(2,0), R[7]=R(2,1),
      //                               R[2]=R(0,2), R[5]=R(1,2)
      const double r33 = R[8];
      double p1 = 0.0;
      double ph = 0.0;
      double p2 = 0.0;

      if(std::abs(std::abs(r33) - 1.0) > k_Eps)
      {
        const double zeta = 1.0 / std::sqrt(1.0 - r33 * r33);
        ph = std::acos(std::clamp(r33, -1.0, 1.0));
        p1 = std::atan2(R[6] * zeta, -R[7] * zeta);
        p2 = std::atan2(R[2] * zeta, R[5] * zeta);
      }
      else if(r33 > 0.0)
      {
        // PHI ≈ 0: degenerate case
        p1 = std::atan2(R[1], R[0]);
        ph = 0.0;
        p2 = 0.0;
      }
      else
      {
        // PHI ≈ π: degenerate case
        p1 = -std::atan2(-R[1], R[0]);
        ph = k_Pi;
        p2 = 0.0;
      }

      // Wrap negative angles into canonical range
      if(p1 < 0.0)
      {
        p1 = std::fmod(p1 + 100.0 * k_TwoPi, k_TwoPi);
      }
      if(ph < 0.0)
      {
        ph = std::fmod(ph + 100.0 * k_Pi, k_Pi);
      }
      if(p2 < 0.0)
      {
        p2 = std::fmod(p2 + 100.0 * k_TwoPi, k_TwoPi);
      }

      const int row = i * numOps + k;
      result(row, 0) = p1;
      result(row, 1) = ph;
      result(row, 2) = p2;
    }
  }

  return result;
}

} // namespace mtrsim

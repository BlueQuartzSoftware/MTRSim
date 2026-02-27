#pragma once

#include <Eigen/Dense>

namespace mtrsim
{

/**
 * @brief Enumerates supported crystal symmetry groups.
 */
enum class CrystalSystem
{
  HCP, ///< Hexagonal close-packed (12 symmetry operators)
  FCC, ///< Face-centred cubic (24 symmetry operators)
};

/**
 * @brief Generates all crystallographically equivalent Euler-angle
 *        representations of a set of orientations.
 *
 * Equivalent to symmetric_euler_angles.m.  Delegates to EbsdLib's
 * LaueOps classes (HexagonalOps / CubicOps) for the symmetry operator
 * tables, using HexagonalOps::getMatSymOpD(i) and
 * OrientationTransformation::eu2om / om2eu for the rotation math.
 */
class CrystalSymmetry
{
public:
  explicit CrystalSymmetry(CrystalSystem system = CrystalSystem::HCP);

  /**
   * @brief Expand orientations to all symmetrically equivalent forms.
   *
   * @param phi1  Input phi1 angles [rad], length N
   * @param phi   Input PHI  angles [rad], length N
   * @param phi2  Input phi2 angles [rad], length N
   * @return      Matrix of shape [N*numOps x 3] with columns [phi1, PHI, phi2]
   */
  Eigen::MatrixXd expand(const Eigen::VectorXd& phi1,
                          const Eigen::VectorXd& phi,
                          const Eigen::VectorXd& phi2) const;

  /// Returns the number of symmetry operators for the configured crystal system.
  int numOperators() const;

private:
  CrystalSystem m_System;
  // Symmetry operators are retrieved on-demand from EbsdLib's LaueOps;
  // no static table is stored here.
};

} // namespace mtrsim

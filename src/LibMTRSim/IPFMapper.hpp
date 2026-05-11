#pragma once

#include "libmtrsim_export.h"
#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

namespace mtrsim {

/**
 * @brief Enumerates the crystal symmetry groups that IPFMapper understands.
 *
 * This is a tiny API-discrimination enum used only to select between the two
 * MATLAB-port colour pipelines (HCP polar mapping vs. cubic). It is NOT a
 * crystallographic abstraction: orientation math always routes through
 * EbsdLib's `LaueOps` and `ebsdlib::CrystalStructure::*` codes directly.
 */
enum class CrystalSystem {
  HCP, ///< Hexagonal close-packed (maps to ebsdlib::CrystalStructure::Hexagonal_High)
  FCC, ///< Face-centred cubic    (maps to ebsdlib::CrystalStructure::Cubic_High)
};

/**
 * @brief Selects the IPF colour-mapping algorithm.
 */
enum class IPFColorScheme {
  EbsdLib, ///< Standard EbsdLib/TSL IPF colouring
  MatLab,  ///< Original MATLAB port colouring (Sparkman 2017)
};

/**
 * @brief RGB colour triple in [0, 255] uint8 range.
 */
struct LIBMTRSIM_EXPORT RGBColor {
  uint8_t r = 0;
  uint8_t g = 0;
  uint8_t b = 0;
};

/**
 * @brief Builds an Inverse Pole Figure (IPF) map from Euler angles.
 *
 * Replaces unit_triangle_IPF_coords.m, IPF_colors.m, and view_IPF_map.m.
 * Two colour-mapping algorithms are available:
 *   - EbsdLib:  delegates to EbsdLib's LaueOps::generateIPFColor().
 *   - MatLab:   faithful port of IPF_colors.m + unit_triangle_IPF_coords.m,
 *               using polar-coordinate colour mapping inside the stereographic
 *               unit triangle.  Currently supports HCP only.
 */
class LIBMTRSIM_EXPORT IPFMapper {
public:
  explicit IPFMapper(CrystalSystem system = CrystalSystem::HCP);

  /**
   * @brief Convert Euler angles directly to IPF RGB colours.
   *
   * @param phi1    phi1 angles [rad], length N
   * @param phi     PHI  angles [rad], length N
   * @param phi2    phi2 angles [rad], length N
   * @param refDir  Reference direction in sample frame (default: Z = {0,0,1}).
   *                Only used by the EbsdLib scheme; the MatLab scheme always
   *                uses the Z-axis.
   * @param scheme  Colour-mapping algorithm (default: EbsdLib)
   * @return        Per-voxel RGB colours, length N, values in [0, 255]
   */
  std::vector<RGBColor>
  eulerToColors(const Eigen::VectorXd &phi1, const Eigen::VectorXd &phi,
                const Eigen::VectorXd &phi2,
                std::array<double, 3> refDir = {0.0, 0.0, 1.0},
                IPFColorScheme scheme = IPFColorScheme::EbsdLib) const;

  /**
   * @brief Render an IPF map image and write it to a PNG file.
   *
   * @param spatialCoords  [N x 2] (x, y) positions [mm]
   * @param phi1           phi1 angles [rad]
   * @param phi            PHI  angles [rad]
   * @param phi2           phi2 angles [rad]
   * @param outputPath     Destination PNG file path
   * @param scheme         Colour-mapping algorithm (default: EbsdLib)
   */
  void writePNG(const Eigen::MatrixXd &spatialCoords,
                const Eigen::VectorXd &phi1, const Eigen::VectorXd &phi,
                const Eigen::VectorXd &phi2, const std::string &outputPath,
                IPFColorScheme scheme = IPFColorScheme::EbsdLib) const;

  /**
   * @brief Render the HCP IPF triangle legend using the MATLAB polar colour
   *        mapping and write it to a PNG file.
   *
   * The image is rectangular, tightly fitting the HCP stereographic unit
   * triangle (the pie-slice from [0001] at the origin to [2-1-10] at (1,0)
   * to [10-10] at (cos 30°, sin 30°), with the circular arc boundary).
   * Pixels outside the triangle are white.
   *
   * Currently supports HCP only (throws for FCC).
   *
   * @param imageDim    Width of the output image in pixels; height is
   *                    computed to preserve the triangle's aspect ratio.
   * @param outputPath  Destination PNG file path.
   */
  void writeIPFTriangleLegendMatLab(int imageDim,
                                    const std::string &outputPath) const;

private:
  CrystalSystem m_System;

  /// EbsdLib colour path (delegates to LaueOps::generateIPFColor).
  std::vector<RGBColor>
  eulerToColorsEbsdLib(const Eigen::VectorXd &phi1, const Eigen::VectorXd &phi,
                       const Eigen::VectorXd &phi2,
                       std::array<double, 3> refDir) const;

  /// MATLAB colour path — port of unit_triangle_IPF_coords.m + IPF_colors.m.
  std::vector<RGBColor> eulerToColorsMatLab(const Eigen::VectorXd &phi1,
                                            const Eigen::VectorXd &phi,
                                            const Eigen::VectorXd &phi2) const;
};

} // namespace mtrsim

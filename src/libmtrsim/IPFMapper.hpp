#pragma once

#include "CrystalSymmetry.hpp"
#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

namespace mtrsim
{

/**
 * @brief RGB colour triple in [0, 255] uint8 range.
 */
struct RGBColor
{
  uint8_t r = 0;
  uint8_t g = 0;
  uint8_t b = 0;
};

/**
 * @brief Builds an Inverse Pole Figure (IPF) map from Euler angles.
 *
 * Replaces unit_triangle_IPF_coords.m, IPF_colors.m, and view_IPF_map.m.
 * Delegates colour computation to EbsdLib's HexagonalOps::generateIPFColor(),
 * which performs the stereographic projection and colour mapping in one step.
 */
class IPFMapper
{
public:
  explicit IPFMapper(CrystalSystem system = CrystalSystem::HCP);

  /**
   * @brief Convert Euler angles directly to IPF RGB colours.
   *
   * Wraps EbsdLib's HexagonalOps::generateIPFColor().  The reference
   * direction is the sample Z-axis (0, 0, 1) by default, matching the
   * MATLAB convention in generate_inverse_pole_figure.m.
   *
   * @param phi1    phi1 angles [rad], length N
   * @param phi     PHI  angles [rad], length N
   * @param phi2    phi2 angles [rad], length N
   * @param refDir  Reference direction in sample frame (default: Z = {0,0,1})
   * @return        Per-voxel RGB colours, length N, values in [0, 255]
   */
  std::vector<RGBColor> eulerToColors(const Eigen::VectorXd& phi1,
                                       const Eigen::VectorXd& phi,
                                       const Eigen::VectorXd& phi2,
                                       std::array<double, 3> refDir = {0.0, 0.0, 1.0}) const;

  /**
   * @brief Render an IPF map image and write it to a PNG file.
   *
   * @param spatialCoords  [N x 2] (x, y) positions [mm]
   * @param phi1           phi1 angles [rad]
   * @param phi            PHI  angles [rad]
   * @param phi2           phi2 angles [rad]
   * @param outputPath     Destination PNG file path
   */
  void writePNG(const Eigen::MatrixXd& spatialCoords,
                const Eigen::VectorXd& phi1,
                const Eigen::VectorXd& phi,
                const Eigen::VectorXd& phi2,
                const std::string& outputPath) const;

private:
  CrystalSystem m_System;
};

} // namespace mtrsim

#pragma once

#include <Eigen/Dense>
#include <string>

namespace mtrsim
{

/**
 * @brief EBSD scan data loaded from CSV files.
 *
 * Matches the variables read by load_MTR_data.m.
 */
struct EBSDData
{
  Eigen::MatrixXd spatialCoords; ///< [N x 2] (X, Y) positions [mm]
  Eigen::MatrixXd eulerAngles;   ///< [N x 3] (phi1, PHI, phi2) [rad]
  Eigen::VectorXi parentIds;     ///< MTR parent grain IDs, length N
  Eigen::VectorX<bool> isMTR;    ///< Boolean MTR membership mask, length N
};

/**
 * @brief Loads experimental EBSD data from a directory containing CSV exports.
 *
 * Expected files in the directory:
 *   - EulerAngles.csv   (N rows × 3 cols: phi1, PHI, phi2 in radians)
 *   - X_Position.csv    (N rows × 1 col)
 *   - Y_Position.csv    (N rows × 1 col)
 *   - ParentIds.csv     (N rows × 1 col, integer grain IDs)
 *   - BoolMTR.csv       (N rows × 1 col, 0 or 1)
 */
class MTRDataLoader
{
public:
  MTRDataLoader() = default;

  /**
   * @brief Load all EBSD CSV files from a directory.
   *
   * @param directoryPath  Path to the directory containing the CSV files
   * @return               Populated EBSDData struct
   */
  EBSDData load(const std::string& directoryPath);
};

} // namespace mtrsim

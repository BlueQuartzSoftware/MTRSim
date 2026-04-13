#pragma once

#include "ODFSampler.hpp"
#include "mtrsim_export.h"
#include <Eigen/Dense>

namespace mtrsim
{

/**
 * @brief Pole figure data: projected (X, Y) coordinates and intensity values.
 */
struct MTRSIM_EXPORT PoleFigureData
{
  Eigen::VectorXd x;         ///< Stereographic X coordinates
  Eigen::VectorXd y;         ///< Stereographic Y coordinates
  Eigen::VectorXd intensity; ///< Normalised intensity per bin
};

/**
 * @brief Converts a discrete ODF to a pole figure.
 *
 * Equivalent to convert_ODF_to_PF.m.  Internally builds an
 * EbsdLib::PoleFigureConfiguration_t from the ODF bins and delegates
 * to HexagonalOps::generatePoleFigure() for the stereographic
 * projection and intensity accumulation.
 */
class MTRSIM_EXPORT PoleFigure
{
public:
  PoleFigure() = default;

  /**
   * @brief Convert an ODF to pole figure intensities.
   *
   * @param component  ODF component (ODFval and bin arrays)
   * @param degSpacing Euler-space bin spacing in degrees (must match ODF)
   * @return           PoleFigureData with stereographic coordinates and intensities
   */
  PoleFigureData fromODF(const ODFComponent& component, double degSpacing = 5.0);
};

} // namespace mtrsim

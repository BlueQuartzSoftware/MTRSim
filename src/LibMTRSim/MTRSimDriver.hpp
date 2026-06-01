#pragma once

#include "libmtrsim_export.h"

#include "ODFSampler.hpp"      // mtrsim::ODFComponent, mtrsim::EulerAngles
#include "SimulationParams.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <random>
#include <vector>

namespace mtrsim {

/**
 * @brief Build a uniform reference ODF on an (n1 x nPHI x n2) Euler grid.
 *
 * The flat bin-centre arrays match exactly what ODFCalculator::compute and the
 * ODFSampler expect, with ix = i1*(nPHI*n2) + iPHI*n2 + i2:
 *   phi1Bins[ix] = (i1   + 0.5) * 2*pi / n1
 *   phiBins[ix]  = (iPHI + 0.5) *   pi / nPHI
 *   phi2Bins[ix] = (i2   + 0.5) * 2*pi / n2
 */
LIBMTRSIM_EXPORT ODFComponent buildUniformODF(int n1, int nPHI, int n2);

/**
 * @brief Reconstruct an ODF component from flat grid data + degree spacing.
 *
 * @param values    Flat ODFval, length n1*nPHI*n2, row-major with
 *                  ix = i1*(nPHI*n2) + iPHI*n2 + i2 (phi1 slowest, phi2 fastest).
 * @param n1,nPHI,n2  Bin counts along phi1, PHI, phi2.
 * @param stepDeg1,stepDegPHI,stepDeg2  Bin sizes [degrees] (geometry spacing).
 * @return ODFComponent with bin centres [rad] and values normalized to sum 1.
 */
LIBMTRSIM_EXPORT ODFComponent gridToODFComponent(const std::vector<double>& values, int n1, int nPHI, int n2, double stepDeg1, double stepDegPHI, double stepDeg2);

} // namespace mtrsim

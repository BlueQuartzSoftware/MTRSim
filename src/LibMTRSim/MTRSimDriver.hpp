#pragma once

#include "libmtrsim_export.h"

#include "ODFSampler.hpp" // mtrsim::ODFComponent, mtrsim::EulerAngles
#include "SimulationParams.hpp"

#include <Eigen/Dense>
#include <cstddef>
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
 *                  ix = i1*(nPHI*n2) + iPHI*n2 + i2 (phi1 slowest, phi2
 * fastest).
 * @param n1,nPHI,n2  Bin counts along phi1, PHI, phi2.
 * @param stepDeg1,stepDegPHI,stepDeg2  Bin sizes [degrees] (geometry spacing).
 * @return ODFComponent with bin centres [rad] and values normalized to sum 1.
 */
LIBMTRSIM_EXPORT ODFComponent
gridToODFComponent(const std::vector<double> &values, int n1, int nPHI, int n2,
                   double stepDeg1, double stepDegPHI, double stepDeg2);

/**
 * @brief Remap a per-voxel vector from simulation order to SIMPLNX z,y,x order.
 *
 * Simulation order:  kSim = iz*(nx*ny) + ix*ny + iy   (y fastest).
 * SIMPLNX order:     kNx  = iz*(ny*nx) + iy*nx + ix   (x fastest).
 *   out[kNx] = in[kSim].
 *
 * @tparam T element type (int or double).
 */
template <typename T>
std::vector<T> remapSimToZYX(const std::vector<T> &in, int nx, int ny, int nz) {
  std::vector<T> out(in.size());
  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {
        const std::size_t kSim = static_cast<std::size_t>(iz) * nx * ny +
                                 static_cast<std::size_t>(ix) * ny + iy;
        const std::size_t kNx = static_cast<std::size_t>(iz) * ny * nx +
                                static_cast<std::size_t>(iy) * nx + ix;
        out[kNx] = in[kSim];
      }
    }
  }
  return out;
}

/**
 * @brief Per-voxel simulation output in SIMPLNX z,y,x order.
 */
struct LIBMTRSIM_EXPORT MTRSimResult {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  std::vector<int32_t> mtrIndex; ///< 1-based component id per voxel, length N
  std::vector<double> phi1;      ///< Euler phi1 [rad] per voxel, length N
  std::vector<double> phi;       ///< Euler PHI  [rad] per voxel, length N
  std::vector<double> phi2;      ///< Euler phi2 [rad] per voxel, length N
};

/**
 * @brief Run the full MTR simulation: PGRF assignment -> per-component ODF
 *        sampling -> per-voxel orientation assignment, returned in SIMPLNX
 *        z,y,x voxel order.
 *
 * @param params         Fully populated SimulationParams (consistent length
 * unit).
 * @param odfComponents  One ODFComponent per volume-fraction entry, shared
 * grid.
 * @param rng            Seeded RNG (mt19937_64).
 * @param n1,nPHI,n2     Bin counts of the ODF grid (for the uniform reference).
 */
LIBMTRSIM_EXPORT MTRSimResult
simulateMTR(const SimulationParams &params,
            const std::vector<ODFComponent> &odfComponents,
            std::mt19937_64 &rng, int n1, int nPHI, int n2);

} // namespace mtrsim

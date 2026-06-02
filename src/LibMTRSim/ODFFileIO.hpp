#pragma once

#include "libmtrsim_export.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace mtrsim
{

/**
 * @brief Metadata describing the contents of a MATLAB-format ODF HDF5 file.
 *
 * The file format stores N ODF components, each as a 3-D regular histogram
 * over Euler space (phi1, PHI, phi2).  All components in a single file share
 * identical bin-edge arrays.
 */
struct LIBMTRSIM_EXPORT ODFFileMetadata
{
  int64_t numComponents;                       ///< Number of ODF components stored in the file (>= 1).
  std::array<int64_t, 3> dimsPhi1PHIPhi2;      ///< Number of bins per axis; phi1
                                               ///< slowest-varying, phi2 fastest.
  std::array<double, 3> spacingDegPhi1PHIPhi2; ///< Uniform bin size in DEGREES, one per axis.
};

/**
 * @brief A single ODF component's raw histogram values.
 *
 * The values buffer is row-major with phi1 slowest-varying and phi2 fastest:
 *   index(i1, iPHI, i2) = i1 * (nPHI * nphi2) + iPHI * nphi2 + i2
 * Size is nphi1 * nPHI * nphi2.
 */
struct LIBMTRSIM_EXPORT ODFFileComponent
{
  std::vector<double> values;
};

/**
 * @brief Read metadata (component count, dims, spacing) without touching ODFval
 * arrays.
 *
 * Validates structural invariants:
 *   - <pathPrefix>/num_components >= 1
 *   - components are present contiguously from component_0 .. component_{N-1}
 *   - every component's phi1_bins / PHI_bins / phi2_bins are byte-exact
 * identical to component_0's arrays
 *   - bin-edge arrays are strictly monotonically increasing and uniformly
 * spaced
 *   - each component's ODFval dataset has the expected size
 *
 * @param file        Path to the HDF5 file.
 * @param pathPrefix  Internal HDF5 group path that contains num_components and
 *                    component_N subgroups. Leading '/' is auto-inserted and
 *                    trailing '/' stripped. Defaults to "/ODF_best" (MATLAB
 *                    convention).
 *
 * @throws std::runtime_error on any validation failure (missing file, missing
 *                            prefix group, missing dataset, inconsistent bins,
 *                            non-uniform spacing, etc.)
 */
LIBMTRSIM_EXPORT ODFFileMetadata readODFMetadata(const std::filesystem::path& file, const std::string& pathPrefix = "/ODF_best");

/**
 * @brief Full read of metadata plus every component's ODFval array.
 *
 * Internally calls readODFMetadata() for full validation before reading
 * any ODFval data.
 *
 * @throws std::runtime_error on any validation or read failure.
 */
LIBMTRSIM_EXPORT std::vector<ODFFileComponent> readODFComponents(const std::filesystem::path& file, const std::string& pathPrefix = "/ODF_best");

/**
 * @brief Write an HDF5 file round-trip-compatible with the MATLAB ODF format.
 *
 * Builds bin-edge arrays in radians from @p dimsPhi1PHIPhi2 and
 * @p spacingDegPhi1PHIPhi2 via edges[i] = i * stepRad for i in [0, N], matching
 * MATLAB's `0:2*pi/num_bins:2*pi` form.
 *
 * @param file                   Output path; overwrites if it exists.
 * @param dimsPhi1PHIPhi2        Number of bins per axis. Each entry must be >
 * 0.
 * @param spacingDegPhi1PHIPhi2  Uniform bin size per axis in DEGREES.
 * @param components             One entry per ODF component; each must have
 *                               values.size() == dims[0] * dims[1] * dims[2].
 * @param pathPrefix             Internal HDF5 group to create and write under.
 *                               Leading '/' is auto-inserted and trailing '/'
 *                               stripped. If the normalized prefix is "/" the
 *                               datasets are written at the file root without
 *                               creating an enclosing group. Defaults to
 *                               "/ODF_best".
 *
 * @throws std::runtime_error on any invalid input or HDF5 error.
 */
LIBMTRSIM_EXPORT void writeODFFile(const std::filesystem::path& file, const std::array<int64_t, 3>& dimsPhi1PHIPhi2, const std::array<double, 3>& spacingDegPhi1PHIPhi2,
                                   const std::vector<ODFFileComponent>& components, const std::string& pathPrefix = "/ODF_best");

/**
 * @brief Read an optional int64 fixture-version field from a reference HDF5.
 *
 * Looks for `<pathPrefix>/fixture_version` in @p file. Returns the value on
 * success, std::nullopt if the field is absent. Tests use this to detect a
 * stale reference (different orientation table) and emit a clear regeneration
 * message instead of producing a confusing bin-by-bin diff.
 *
 * @param file         Path to the HDF5 file.
 * @param pathPrefix   Internal HDF5 group path. Defaults to "/ODF_best".
 *
 * @throws std::runtime_error if the file can't be opened or the field exists
 *         but isn't an int64 scalar.
 */
LIBMTRSIM_EXPORT std::optional<int64_t> tryReadFixtureVersion(const std::filesystem::path& file, const std::string& pathPrefix = "/ODF_best");

} // namespace mtrsim

#include "ODFFileIO.hpp"

#include <hdf5.h>

#include <cmath>
#include <cstring>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mtrsim {

namespace {

// ─────────────────────────────────────────────────────────────────────────────
// RAII closers for HDF5 handles. Each holds an hid_t, closes on destruction
// if still open, and exposes release()/reset() for explicit hand-off.
// ─────────────────────────────────────────────────────────────────────────────

template <herr_t (*CloseFn)(hid_t)> class H5Handle {
public:
  H5Handle() = default;

  explicit H5Handle(hid_t id) : m_Id(id) {}

  H5Handle(const H5Handle &) = delete;
  H5Handle &operator=(const H5Handle &) = delete;

  H5Handle(H5Handle &&other) noexcept : m_Id(other.m_Id) { other.m_Id = -1; }

  H5Handle &operator=(H5Handle &&other) noexcept {
    if (this != &other) {
      reset();
      m_Id = other.m_Id;
      other.m_Id = -1;
    }
    return *this;
  }

  ~H5Handle() { reset(); }

  void reset(hid_t newId = -1) {
    if (m_Id >= 0) {
      CloseFn(m_Id);
    }
    m_Id = newId;
  }

  hid_t get() const { return m_Id; }

  bool valid() const { return m_Id >= 0; }

private:
  hid_t m_Id{-1};
};

using FileHandle = H5Handle<&H5Fclose>;
using GroupHandle = H5Handle<&H5Gclose>;
using DatasetHandle = H5Handle<&H5Dclose>;
using DataspaceHandle = H5Handle<&H5Sclose>;

// ─────────────────────────────────────────────────────────────────────────────
// RAII scope guard: save the current HDF5 auto-error-printing state, disable
// it for the duration of the scope, and restore the previous callback in the
// destructor. We surface failures via exceptions with our own messages, so we
// don't want the library smearing its own error stack to stderr in the middle
// of an expected failure path (e.g. probing a group that we know might not
// exist).
// ─────────────────────────────────────────────────────────────────────────────
class H5ErrorSuppressor {
public:
  H5ErrorSuppressor() {
    H5Eget_auto2(H5E_DEFAULT, &m_OldFunc, &m_OldClientData);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  }
  ~H5ErrorSuppressor() {
    H5Eset_auto2(H5E_DEFAULT, m_OldFunc, m_OldClientData);
  }
  H5ErrorSuppressor(const H5ErrorSuppressor &) = delete;
  H5ErrorSuppressor(H5ErrorSuppressor &&) = delete;
  H5ErrorSuppressor &operator=(const H5ErrorSuppressor &) = delete;
  H5ErrorSuppressor &operator=(H5ErrorSuppressor &&) = delete;

private:
  H5E_auto2_t m_OldFunc = nullptr;
  void *m_OldClientData = nullptr;
};

// ─────────────────────────────────────────────────────────────────────────────
// Low-level HDF5 read/write helpers
// ─────────────────────────────────────────────────────────────────────────────

int64_t readScalarInt64(hid_t loc, const std::string &path) {
  DatasetHandle ds{H5Dopen2(loc, path.c_str(), H5P_DEFAULT)};
  if (!ds.valid()) {
    throw std::runtime_error("HDF5: cannot open dataset: " + path);
  }
  int64_t value = 0;
  if (H5Dread(ds.get(), H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              &value) < 0) {
    throw std::runtime_error("HDF5: failed to read scalar int64 dataset: " +
                             path);
  }
  return value;
}

std::vector<double> readDoubleVector(hid_t loc, const std::string &path) {
  DatasetHandle ds{H5Dopen2(loc, path.c_str(), H5P_DEFAULT)};
  if (!ds.valid()) {
    throw std::runtime_error("HDF5: cannot open dataset: " + path);
  }

  DataspaceHandle space{H5Dget_space(ds.get())};
  if (!space.valid()) {
    throw std::runtime_error("HDF5: cannot get dataspace for dataset: " + path);
  }

  const int ndims = H5Sget_simple_extent_ndims(space.get());
  if (ndims < 0) {
    throw std::runtime_error("HDF5: failed to get rank of dataset: " + path);
  }
  if (ndims < 1) {
    throw std::runtime_error("HDF5: dataset has invalid rank: " + path);
  }
  std::vector<hsize_t> dims(static_cast<std::size_t>(ndims));
  if (H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr) < 0) {
    throw std::runtime_error("HDF5: failed to read dims of dataset: " + path);
  }

  hsize_t total = 1;
  for (hsize_t d : dims) {
    total *= d;
  }

  std::vector<double> buf(total);
  if (H5Dread(ds.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              buf.data()) < 0) {
    throw std::runtime_error("HDF5: failed to read double dataset: " + path);
  }
  return buf;
}

void writeScalarInt64(hid_t loc, const std::string &name, int64_t value) {
  DataspaceHandle space{H5Screate(H5S_SCALAR)};
  if (!space.valid()) {
    throw std::runtime_error("HDF5: failed to create scalar dataspace for: " +
                             name);
  }
  DatasetHandle ds{H5Dcreate2(loc, name.c_str(), H5T_NATIVE_INT64, space.get(),
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
  if (!ds.valid()) {
    throw std::runtime_error("HDF5: failed to create scalar int64 dataset: " +
                             name);
  }
  if (H5Dwrite(ds.get(), H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               &value) < 0) {
    throw std::runtime_error("HDF5: failed to write scalar int64 dataset: " +
                             name);
  }
}

void writeDoubleVector(hid_t loc, const std::string &name,
                       const std::vector<double> &buf) {
  const hsize_t dims[1] = {static_cast<hsize_t>(buf.size())};
  DataspaceHandle space{H5Screate_simple(1, dims, nullptr)};
  if (!space.valid()) {
    throw std::runtime_error("HDF5: failed to create simple dataspace for: " +
                             name);
  }
  DatasetHandle ds{H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, space.get(),
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
  if (!ds.valid()) {
    throw std::runtime_error("HDF5: failed to create double dataset: " + name);
  }
  if (H5Dwrite(ds.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               buf.data()) < 0) {
    throw std::runtime_error("HDF5: failed to write double dataset: " + name);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation helpers
// ─────────────────────────────────────────────────────────────────────────────

constexpr double k_SpacingTolRad = 1.0e-12;

// Verify a bin-edge array is strictly monotonically increasing and uniformly
// spaced (step consistent within k_SpacingTolRad). Returns the uniform step in
// radians. Throws on any violation.
double validateUniformBins(const std::vector<double> &edges,
                           const std::string &axisName) {
  if (edges.size() < 2) {
    throw std::runtime_error("ODFFileIO: bin-edge array too short on axis: " +
                             axisName);
  }
  const double step0 = edges[1] - edges[0];
  if (step0 <= 0.0) {
    throw std::runtime_error(
        "ODFFileIO: bin edges not strictly increasing on axis: " + axisName);
  }
  for (std::size_t i = 1; i + 1 < edges.size(); ++i) {
    const double step = edges[i + 1] - edges[i];
    if (step <= 0.0) {
      throw std::runtime_error(
          "ODFFileIO: bin edges not strictly increasing on axis: " + axisName);
    }
    if (std::abs(step - step0) > k_SpacingTolRad) {
      std::ostringstream oss;
      oss << "ODFFileIO: non-uniform bin spacing on axis " << axisName
          << " at index " << i << " (step=" << step << ", expected=" << step0
          << ", tol=" << k_SpacingTolRad << ")";
      throw std::runtime_error(oss.str());
    }
  }
  return step0;
}

// Byte-exact compare of two double arrays. Sizes must match first.
bool byteEqual(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size()) {
    return false;
  }
  return std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

// Normalize a user-supplied HDF5 path prefix:
//   - Insert a leading '/' if missing.
//   - Strip a trailing '/' unless the result would be empty.
//   - Empty input and "/" both collapse to "/".
// Examples:
//   ""            -> "/"
//   "/"           -> "/"
//   "ODF_best"    -> "/ODF_best"
//   "/ODF_best/"  -> "/ODF_best"
//   "/a/b/"       -> "/a/b"
std::string normalizePathPrefix(const std::string &prefix) {
  if (prefix.empty()) {
    return "/";
  }
  std::string out = prefix;
  if (out.front() != '/') {
    out.insert(out.begin(), '/');
  }
  while (out.size() > 1 && out.back() == '/') {
    out.pop_back();
  }
  return out;
}

std::string componentGroupName(const std::string &normalizedPrefix,
                               int64_t idx) {
  if (normalizedPrefix == "/") {
    return "/component_" + std::to_string(idx);
  }
  return normalizedPrefix + "/component_" + std::to_string(idx);
}

std::string joinPath(const std::string &normalizedPrefix,
                     const std::string &leaf) {
  if (normalizedPrefix == "/") {
    return "/" + leaf;
  }
  return normalizedPrefix + "/" + leaf;
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

ODFFileMetadata readODFMetadata(const std::filesystem::path &file,
                                const std::string &pathPrefix) {
  if (!std::filesystem::exists(file)) {
    throw std::runtime_error("ODFFileIO: file does not exist: " +
                             file.string());
  }

  const std::string prefix = normalizePathPrefix(pathPrefix);

  // Suppress HDF5's default stderr error printing; we surface errors via
  // exceptions with our own messages.
  H5ErrorSuppressor suppressErrors;

  FileHandle f{H5Fopen(file.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
  if (!f.valid()) {
    throw std::runtime_error("ODFFileIO: cannot open HDF5 file: " +
                             file.string());
  }

  // Verify the prefix group exists (skip this check when prefix == "/" which is
  // always valid — it's the implicit root).
  if (prefix != "/") {
    const htri_t prefixExists = H5Lexists(f.get(), prefix.c_str(), H5P_DEFAULT);
    if (prefixExists < 0) {
      throw std::runtime_error("ODFFileIO: HDF5 error probing '" + prefix +
                               "' in " + file.string());
    }
    if (prefixExists == 0) {
      throw std::runtime_error("ODFFileIO: group '" + prefix +
                               "' not found in file " + file.string());
    }
  }

  const int64_t numComponents =
      readScalarInt64(f.get(), joinPath(prefix, "num_components"));
  if (numComponents < 1) {
    throw std::runtime_error("ODFFileIO: num_components must be >= 1, got " +
                             std::to_string(numComponents));
  }

  // Verify every expected component group exists
  for (int64_t j = 0; j < numComponents; ++j) {
    const std::string group = componentGroupName(prefix, j);
    const htri_t exists = H5Lexists(f.get(), group.c_str(), H5P_DEFAULT);
    if (exists < 0) {
      throw std::runtime_error("ODFFileIO: HDF5 error probing '" + group +
                               "' in " + file.string());
    }
    if (exists == 0) {
      throw std::runtime_error("ODFFileIO: missing component group: " + group);
    }
  }

  // Read component_0's bin arrays as the reference
  const std::string c0 = componentGroupName(prefix, 0);
  std::vector<double> refPhi1 = readDoubleVector(f.get(), c0 + "/phi1_bins");
  std::vector<double> refPHI = readDoubleVector(f.get(), c0 + "/PHI_bins");
  std::vector<double> refPhi2 = readDoubleVector(f.get(), c0 + "/phi2_bins");

  const double stepPhi1 = validateUniformBins(refPhi1, "phi1");
  const double stepPHI = validateUniformBins(refPHI, "PHI");
  const double stepPhi2 = validateUniformBins(refPhi2, "phi2");

  const int64_t nphi1 = static_cast<int64_t>(refPhi1.size()) - 1;
  const int64_t nPHI = static_cast<int64_t>(refPHI.size()) - 1;
  const int64_t nphi2 = static_cast<int64_t>(refPhi2.size()) - 1;
  const std::size_t expectedOdfSize = static_cast<std::size_t>(nphi1) *
                                      static_cast<std::size_t>(nPHI) *
                                      static_cast<std::size_t>(nphi2);

  // Validate all components: bins byte-identical, ODFval size correct
  for (int64_t j = 0; j < numComponents; ++j) {
    const std::string base = componentGroupName(prefix, j);

    if (j > 0) {
      const std::vector<double> phi1Arr =
          readDoubleVector(f.get(), base + "/phi1_bins");
      const std::vector<double> PHIArr =
          readDoubleVector(f.get(), base + "/PHI_bins");
      const std::vector<double> phi2Arr =
          readDoubleVector(f.get(), base + "/phi2_bins");

      if (!byteEqual(phi1Arr, refPhi1)) {
        throw std::runtime_error(
            "ODFFileIO: phi1_bins differs from component_0 in " + base);
      }
      if (!byteEqual(PHIArr, refPHI)) {
        throw std::runtime_error(
            "ODFFileIO: PHI_bins differs from component_0 in " + base);
      }
      if (!byteEqual(phi2Arr, refPhi2)) {
        throw std::runtime_error(
            "ODFFileIO: phi2_bins differs from component_0 in " + base);
      }
    }

    const std::string odfvalPath = base + "/ODFval";
    DatasetHandle ds{H5Dopen2(f.get(), odfvalPath.c_str(), H5P_DEFAULT)};
    if (!ds.valid()) {
      throw std::runtime_error("ODFFileIO: missing ODFval in " + base);
    }
    DataspaceHandle space{H5Dget_space(ds.get())};
    if (!space.valid()) {
      throw std::runtime_error(
          "ODFFileIO: cannot get dataspace for ODFval in " + base);
    }
    const int ndims = H5Sget_simple_extent_ndims(space.get());
    if (ndims < 0) {
      throw std::runtime_error("ODFFileIO: failed to get rank of dataset '" +
                               odfvalPath + "' in " + file.string());
    }
    std::vector<hsize_t> dims(static_cast<std::size_t>(ndims));
    if (H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr) < 0) {
      throw std::runtime_error("ODFFileIO: failed to read dims of dataset '" +
                               odfvalPath + "' in " + file.string());
    }
    hsize_t total = 1;
    for (hsize_t d : dims) {
      total *= d;
    }
    if (total != expectedOdfSize) {
      std::ostringstream oss;
      oss << "ODFFileIO: ODFval size mismatch in " << base << " (got " << total
          << ", expected " << expectedOdfSize << ")";
      throw std::runtime_error(oss.str());
    }
  }

  constexpr double k_RadToDeg = 180.0 / std::numbers::pi;

  ODFFileMetadata md;
  md.numComponents = numComponents;
  md.dimsPhi1PHIPhi2 = {nphi1, nPHI, nphi2};
  md.spacingDegPhi1PHIPhi2 = {stepPhi1 * k_RadToDeg, stepPHI * k_RadToDeg,
                              stepPhi2 * k_RadToDeg};
  return md;
}

std::vector<ODFFileComponent>
readODFComponents(const std::filesystem::path &file,
                  const std::string &pathPrefix) {
  const std::string prefix = normalizePathPrefix(pathPrefix);

  // Full validation up front — throws if anything is wrong.
  const ODFFileMetadata md = readODFMetadata(file, prefix);

  // Re-suppress HDF5 default error printing for our own file ops
  H5ErrorSuppressor suppressErrors;

  FileHandle f{H5Fopen(file.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
  if (!f.valid()) {
    throw std::runtime_error("ODFFileIO: cannot open HDF5 file: " +
                             file.string());
  }

  const std::size_t expectedSize =
      static_cast<std::size_t>(md.dimsPhi1PHIPhi2[0]) *
      static_cast<std::size_t>(md.dimsPhi1PHIPhi2[1]) *
      static_cast<std::size_t>(md.dimsPhi1PHIPhi2[2]);

  std::vector<ODFFileComponent> out;
  out.reserve(static_cast<std::size_t>(md.numComponents));
  for (int64_t j = 0; j < md.numComponents; ++j) {
    const std::string base = componentGroupName(prefix, j);
    std::vector<double> buf = readDoubleVector(f.get(), base + "/ODFval");
    if (buf.size() != expectedSize) {
      // readODFMetadata should already have caught this, but verify.
      std::ostringstream oss;
      oss << "ODFFileIO: ODFval size mismatch in " << base << " (got "
          << buf.size() << ", expected " << expectedSize << ")";
      throw std::runtime_error(oss.str());
    }
    ODFFileComponent comp;
    comp.values = std::move(buf);
    out.push_back(std::move(comp));
  }

  return out;
}

void writeODFFile(const std::filesystem::path &file,
                  const std::array<int64_t, 3> &dimsPhi1PHIPhi2,
                  const std::array<double, 3> &spacingDegPhi1PHIPhi2,
                  const std::vector<ODFFileComponent> &components,
                  const std::string &pathPrefix) {
  if (components.empty()) {
    throw std::runtime_error(
        "ODFFileIO::writeODFFile: components vector is empty");
  }

  const std::string prefix = normalizePathPrefix(pathPrefix);
  for (std::size_t i = 0; i < 3; ++i) {
    if (dimsPhi1PHIPhi2[i] <= 0) {
      throw std::runtime_error("ODFFileIO::writeODFFile: dims must all be > 0");
    }
  }
  const std::size_t expectedSize =
      static_cast<std::size_t>(dimsPhi1PHIPhi2[0]) *
      static_cast<std::size_t>(dimsPhi1PHIPhi2[1]) *
      static_cast<std::size_t>(dimsPhi1PHIPhi2[2]);
  for (std::size_t ci = 0; ci < components.size(); ++ci) {
    if (components[ci].values.size() != expectedSize) {
      std::ostringstream oss;
      oss << "ODFFileIO::writeODFFile: component " << ci
          << " has values.size()=" << components[ci].values.size()
          << ", expected " << expectedSize;
      throw std::runtime_error(oss.str());
    }
  }

  // Suppress default HDF5 error printing
  H5ErrorSuppressor suppressErrors;

  // Build edge arrays (in radians) from dims + deg spacing:
  //   edges[i] = i * stepRad  for i in [0, N]
  constexpr double k_DegToRad = std::numbers::pi / 180.0;
  auto buildEdges = [](int64_t n, double stepDeg) {
    std::vector<double> edges(static_cast<std::size_t>(n) + 1);
    const double stepRad = stepDeg * k_DegToRad;
    for (std::size_t i = 0; i < edges.size(); ++i) {
      edges[i] = static_cast<double>(i) * stepRad;
    }
    return edges;
  };
  const std::vector<double> phi1Edges =
      buildEdges(dimsPhi1PHIPhi2[0], spacingDegPhi1PHIPhi2[0]);
  const std::vector<double> PHIEdges =
      buildEdges(dimsPhi1PHIPhi2[1], spacingDegPhi1PHIPhi2[1]);
  const std::vector<double> phi2Edges =
      buildEdges(dimsPhi1PHIPhi2[2], spacingDegPhi1PHIPhi2[2]);

  FileHandle f{H5Fcreate(file.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT)};
  if (!f.valid()) {
    throw std::runtime_error(
        "ODFFileIO::writeODFFile: cannot create HDF5 file: " + file.string());
  }

  // When the normalized prefix is "/" we write directly to the file root and
  // skip the enclosing group. Otherwise we create the prefix group once and
  // write num_components + component_N subgroups inside it.
  GroupHandle odfGroup;
  hid_t containerId = f.get();
  if (prefix != "/") {
    odfGroup.reset(H5Gcreate2(f.get(), prefix.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                              H5P_DEFAULT));
    if (!odfGroup.valid()) {
      throw std::runtime_error("ODFFileIO::writeODFFile: cannot create group " +
                               prefix);
    }
    containerId = odfGroup.get();
  }

  writeScalarInt64(containerId, "num_components",
                   static_cast<int64_t>(components.size()));

  for (std::size_t j = 0; j < components.size(); ++j) {
    const std::string name = "component_" + std::to_string(j);
    GroupHandle compGroup{H5Gcreate2(containerId, name.c_str(), H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT)};
    if (!compGroup.valid()) {
      throw std::runtime_error("ODFFileIO::writeODFFile: cannot create group " +
                               name);
    }
    writeDoubleVector(compGroup.get(), "ODFval", components[j].values);
    writeDoubleVector(compGroup.get(), "phi1_bins", phi1Edges);
    writeDoubleVector(compGroup.get(), "PHI_bins", PHIEdges);
    writeDoubleVector(compGroup.get(), "phi2_bins", phi2Edges);
  }
}

std::optional<int64_t> tryReadFixtureVersion(const std::filesystem::path &file,
                                             const std::string &pathPrefix) {
  if (!std::filesystem::exists(file)) {
    throw std::runtime_error(
        "ODFFileIO::tryReadFixtureVersion: file does not exist: " +
        file.string());
  }

  const std::string prefix = normalizePathPrefix(pathPrefix);
  const std::string versionPath = joinPath(prefix, "fixture_version");

  H5ErrorSuppressor suppressErrors;

  FileHandle f{H5Fopen(file.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
  if (!f.valid()) {
    throw std::runtime_error(
        "ODFFileIO::tryReadFixtureVersion: cannot open HDF5 file: " +
        file.string());
  }

  // Field is optional: probe with H5Lexists and return nullopt if absent.
  const htri_t exists = H5Lexists(f.get(), versionPath.c_str(), H5P_DEFAULT);
  if (exists < 0) {
    throw std::runtime_error(
        "ODFFileIO::tryReadFixtureVersion: HDF5 error probing '" + versionPath +
        "' in " + file.string());
  }
  if (exists == 0) {
    return std::nullopt;
  }
  return readScalarInt64(f.get(), versionPath);
}

} // namespace mtrsim

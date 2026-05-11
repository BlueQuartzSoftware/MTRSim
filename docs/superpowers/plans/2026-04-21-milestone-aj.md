# Milestone AJ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver the SIMPLNX plugin integration for MTRSim — three filters (Import, Export, ComputeODF), two new `LibMTRSim` library helpers (`ODFBuilder`, `SymmetricEulers`), a `LibMTRSim`-internal HDF5 I/O helper, a dual-build plugin/library CMake shape, and the full test suite.

**Architecture:** A new `MTRSim` simplnx plugin is hosted in the existing MTRSim repository. The plugin's default build compiles `LibMTRSim` sources directly into the plugin target. A `MTRSIM_BUILD_STANDALONE_LIB` option opts-in to additionally building `LibMTRSim` as a shared library, the `mtrsim` CLI, and the Python bindings — preserving the current MATLAB-comparison workflow. All three filters share the "ODF as ImageGeom" DataStructure convention (see spec Section 3).

**Tech Stack:** C++20, CMake 3.26+, vcpkg, simplnx core, EbsdLib, HDF5 1.14+, Catch2 v3 (tests), Python 3 (scaffolding via `make_filter.py`).

**Design spec:** `docs/superpowers/specs/2026-04-21-milestone-aj-design.md`

**Dependency graph:**
```
T1 (scaffolding) ─┐
T2 (ODFFileIO) ───┴─→ T3 (Import) ──→ T4 (Export + round-trip)
T5 (SymmetricEulers) ─→ T6 (ODFBuilder) ─→ T7 (ComputeODF filter) ─→ T8 (docs)
```
T2 and T5 are independent of everything and of each other — good subagent-parallelism targets.

---

## File Structure

**New files:**

```
CMakeLists.txt                                  (rewrite for dual-build)
src/MTRSim/CMakeLists.txt                       (plugin subdir CMakeLists)
src/MTRSim/src/MTRSim/MTRSimPlugin.hpp
src/MTRSim/src/MTRSim/MTRSimPlugin.cpp
src/MTRSim/src/MTRSim/Filters/
    ReadMTRSimODFFilter.{hpp,cpp}
    WriteMTRSimODFFilter.{hpp,cpp}
    ComputeODFFilter.{hpp,cpp}
    Algorithms/
      ReadMTRSimODF.{hpp,cpp}
      WriteMTRSimODF.{hpp,cpp}
      ComputeODF.{hpp,cpp}
src/LibMTRSim/ODFFileIO.{hpp,cpp}               (HDF5 read/write for MATLAB ODF format)
src/LibMTRSim/SymmetricEulers.{hpp,cpp}         (EbsdLib symmetry expansion)
src/LibMTRSim/ODFBuilder.{hpp,cpp}              (binning + smoothing accumulator)
tests/test_odf_file_io.cpp
tests/test_symmetric_eulers.cpp
tests/test_odf_builder.cpp
test/ReadMTRSimODFFilterTest.cpp
test/WriteMTRSimODFFilterTest.cpp
test/ComputeODFFilterTest.cpp
test/test_data/
    simulation_ODF_roundtrip.h5.gz              (small exemplar from data/simulation_ODF.h5)
    hcp_euler_sample.h5.gz                      (small EBSD fixture)
    hcp_euler_sample_odf_reference.h5.gz        (MATLAB-computed reference ODF)
docs/ReadMTRSimODFFilter.md
docs/WriteMTRSimODFFilter.md
docs/ComputeODFFilter.md
```

**Modified files:**

```
src/LibMTRSim/CMakeLists.txt   (add new helper sources)
tests/CMakeLists.txt            (register new library-level tests)
vcpkg.json                      (no new deps expected; verify Catch2 + HDF5)
```

---

## Task 1: Plugin scaffolding & dual-build CMake

**Files:**
- Modify: `CMakeLists.txt`
- Create: `src/MTRSim/CMakeLists.txt`
- Create: `src/MTRSim/src/MTRSim/MTRSimPlugin.{hpp,cpp}`
- Run: `make_filter.py` three times (generates filter skeletons; we fill them in later tasks)

### Context

The existing top-level `CMakeLists.txt` builds `LibMTRSim` (`add_subdirectory(src/LibMTRSim)`), `app/`, `tools/`, `tests/`, and `wrapping/python/`. We're restructuring so that the simplnx `MTRSim` plugin is the default build target, while the library/CLI/Python stack becomes opt-in.

simplnx plugin CMake convention: each plugin includes `${simplnx_SOURCE_DIR}/cmake/Plugin.cmake`, declares a `FilterList`, and invokes `create_simplnx_plugin(...)`. Reference example: `/Users/mjackson/Workspace7/simplnx/src/Plugins/SimplnxCore/CMakeLists.txt`. This plan does NOT enumerate the CMake details — the project lead (Mike Jackson) authors them. This task establishes **the target shape and verification gates**.

- [ ] **Step 1: Add `MTRSIM_BUILD_STANDALONE_LIB` option to top-level CMakeLists.txt**

Insert after the existing option blocks (around line 60):

```cmake
# ------------------------------------------------------------------------------
# Dual-build mode: when ON, build standalone LibMTRSim + CLI + Python bindings
# in addition to the simplnx plugin. Default OFF — plugin-only build.
# ------------------------------------------------------------------------------
option(MTRSIM_BUILD_STANDALONE_LIB "Build the standalone LibMTRSim shared library, the mtrsim CLI, and Python bindings" OFF)
```

Gate the existing `add_subdirectory(app)`, `add_subdirectory(tools)`, and `add_subdirectory(wrapping/python)` blocks behind `if(MTRSIM_BUILD_STANDALONE_LIB)`. Keep `add_subdirectory(src/LibMTRSim)` unconditional — the plugin will link against it.

- [ ] **Step 2: Add plugin subdirectory**

At the end of `CMakeLists.txt` (after the existing sub-projects block):

```cmake
# ------------------------------------------------------------------------------
# SIMPLNX plugin — always built. Requires simplnx as a parent project or via
# `SIMPLNX_SOURCE_DIR` being set. Plugin target depends on LibMTRSim.
# ------------------------------------------------------------------------------
add_subdirectory(src/MTRSim)
```

- [ ] **Step 3: Create plugin-level CMakeLists.txt**

Follow the `SimplnxCore` pattern:

```cmake
# src/MTRSim/CMakeLists.txt
include("${simplnx_SOURCE_DIR}/cmake/Plugin.cmake")

set(PLUGIN_NAME "MTRSim")
set(${PLUGIN_NAME}_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})

set(FilterList
  ReadMTRSimODFFilter
  WriteMTRSimODFFilter
  ComputeODFFilter
)

set(ActionList
  # none for AJ
)

set(AlgorithmList
  ReadMTRSimODF
  WriteMTRSimODF
  ComputeODF
)

create_simplnx_plugin(NAME ${PLUGIN_NAME}
                     FILTER_LIST ${FilterList}
                     ACTION_LIST ${ActionList}
                     ALGORITHM_LIST ${AlgorithmList}
                     DESCRIPTION "MTRSim - MTR Virtual Representation Generation"
                     VERSION "1.0.0"
                     DOC_CHECK
)

target_link_libraries(${PLUGIN_NAME} PRIVATE mtrsim)
```

(The `simplnx` macro `create_simplnx_plugin` generates the plugin target, headers, registration code, and test CMakeLists. Exact macro signature may differ — reference `SimplnxCore/CMakeLists.txt` and adjust.)

- [ ] **Step 4: Generate filter skeletons via `make_filter.py`**

Run these three commands:

```bash
cd /Users/mjackson/Workspace7/MTRSim
SIMPLNX=/Users/mjackson/Workspace7/simplnx
mkdir -p src/MTRSim/src/MTRSim/Filters/Algorithms
mkdir -p docs
mkdir -p src/MTRSim/test

python "$SIMPLNX/scripts/make_filter.py" \
  --plugin_dir "$PWD/src/MTRSim" \
  --name ReadMTRSimODF \
  --template_dir "$SIMPLNX/scripts"

python "$SIMPLNX/scripts/make_filter.py" \
  --plugin_dir "$PWD/src/MTRSim" \
  --name WriteMTRSimODF \
  --template_dir "$SIMPLNX/scripts"

python "$SIMPLNX/scripts/make_filter.py" \
  --plugin_dir "$PWD/src/MTRSim" \
  --name ComputeODF \
  --template_dir "$SIMPLNX/scripts"
```

Each invocation writes `{Name}Filter.{hpp,cpp}`, `Algorithms/{Name}.{hpp,cpp}`, a test stub, and a docs stub. Inspect them — subsequent tasks edit the TODO markers.

- [ ] **Step 5: Verify configure**

From a scratch build directory, configure with simplnx + MTRSim and verify the plugin is listed:

```bash
cd /Users/mjackson/Workspace7/Build
cmake --preset mtrsim-Rel 2>&1 | tee configure.log
grep -E "MTRSim|Plugin" configure.log | head -20
```

**Expected:** "Enabling plugin: MTRSim" (or similar simplnx log line). No CMake errors.

- [ ] **Step 6: Verify plugin builds (empty filter bodies)**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target MTRSim 2>&1 | tail -20
```

**Expected:** Plugin shared library built, linked against `mtrsim` (LibMTRSim). No symbol errors. The three filter skeletons compile as no-op stubs because `make_filter.py` generates compilable boilerplate.

- [ ] **Step 7: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add CMakeLists.txt src/MTRSim
git commit -m "$(cat <<'EOF'
feat(plugin): scaffold MTRSim simplnx plugin with dual-build CMake

Add MTRSIM_BUILD_STANDALONE_LIB option (default OFF) to gate the
existing lib+CLI+Python stack. Plugin is the default build target; it
links LibMTRSim sources directly. Three empty filter skeletons
generated via make_filter.py.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `mtrsim::ODFFileIO` library helper (TDD)

**Files:**
- Create: `src/LibMTRSim/ODFFileIO.hpp`
- Create: `src/LibMTRSim/ODFFileIO.cpp`
- Create: `tests/test_odf_file_io.cpp`
- Modify: `src/LibMTRSim/CMakeLists.txt` (add new sources)
- Modify: `tests/CMakeLists.txt` (register test)

### Context

HDF5 read/write helpers for the MATLAB-compatible ODF format. Implements spec Section 4 (`readODFMetadata`, `readODFComponents`) plus a writer. The layout is documented in `app/main.cpp:60-65`. No simplnx dependency here — pure `LibMTRSim` + HDF5. This is the best candidate for parallel subagent execution alongside Task 5.

- [ ] **Step 1: Write the header interface**

Create `src/LibMTRSim/ODFFileIO.hpp`:

```cpp
#pragma once

#include "libmtrsim_export.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mtrsim {

/// Metadata for a MATLAB-compatible ODF HDF5 file.
/// Dims are in (phi1, PHI, phi2) order (MATLAB / row-major-flat convention).
/// spacingDeg is in the same order.
struct MTRSIM_EXPORT ODFFileMetadata {
  int64_t                numComponents;
  std::array<int64_t, 3> dimsPhi1PHIPhi2;       // bin counts: phi1 (slowest), PHI, phi2 (fastest)
  std::array<double, 3>  spacingDegPhi1PHIPhi2; // bin size in degrees
};

/// A single component's ODFval payload.
struct MTRSIM_EXPORT ODFFileComponent {
  std::vector<double> values;  // row-major, size = nphi1 * nPHI * nphi2
};

/// Reads /ODF_best/num_components and per-component bin arrays.
/// Validates that every component has identical bin arrays (byte-exact).
/// Throws std::runtime_error with a descriptive message on failure.
MTRSIM_EXPORT ODFFileMetadata readODFMetadata(const std::filesystem::path& file);

/// Reads the full file including all ODFval arrays. Returns one component
/// payload per HDF5 group `component_0..component_{N-1}`.
MTRSIM_EXPORT std::vector<ODFFileComponent> readODFComponents(const std::filesystem::path& file);

/// Writes a MATLAB-compatible ODF HDF5 file. `components` holds N ODFval
/// payloads in row-major order; `dims` and `spacingDeg` define the bin grid.
/// The phi1_bins / PHI_bins / phi2_bins edge arrays are generated from the
/// spacing + dims and converted deg->rad before writing.
MTRSIM_EXPORT void writeODFFile(const std::filesystem::path& file,
                                const std::array<int64_t, 3>& dimsPhi1PHIPhi2,
                                const std::array<double, 3>& spacingDegPhi1PHIPhi2,
                                const std::vector<ODFFileComponent>& components);

}  // namespace mtrsim
```

- [ ] **Step 2: Write the failing test for `readODFMetadata`**

Create `tests/test_odf_file_io.cpp`:

```cpp
#include "LibMTRSim/ODFFileIO.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <filesystem>

using Catch::Approx;
namespace fs = std::filesystem;

// Fixture path relative to the test binary working directory.
// CMake target_compile_definitions sets MTRSIM_TEST_DATA_DIR.
static fs::path testDataDir()
{
  return fs::path(MTRSIM_TEST_DATA_DIR);
}

TEST_CASE("ODFFileIO::readODFMetadata reports correct counts and spacing", "[ODFFileIO]")
{
  const fs::path odf = testDataDir() / "simulation_ODF.h5";
  REQUIRE(fs::exists(odf));

  const auto meta = mtrsim::readODFMetadata(odf);

  REQUIRE(meta.numComponents == 3);
  REQUIRE(meta.dimsPhi1PHIPhi2[0] == 72);  // phi1
  REQUIRE(meta.dimsPhi1PHIPhi2[1] == 36);  // PHI
  REQUIRE(meta.dimsPhi1PHIPhi2[2] == 72);  // phi2
  REQUIRE(meta.spacingDegPhi1PHIPhi2[0] == Approx(5.0).epsilon(1e-9));
  REQUIRE(meta.spacingDegPhi1PHIPhi2[1] == Approx(5.0).epsilon(1e-9));
  REQUIRE(meta.spacingDegPhi1PHIPhi2[2] == Approx(5.0).epsilon(1e-9));
}

TEST_CASE("ODFFileIO::readODFMetadata rejects missing file", "[ODFFileIO]")
{
  REQUIRE_THROWS_AS(mtrsim::readODFMetadata("/nonexistent.h5"), std::runtime_error);
}
```

- [ ] **Step 3: Register the test in tests/CMakeLists.txt**

Append to `tests/CMakeLists.txt`:

```cmake
add_executable(test_odf_file_io test_odf_file_io.cpp)
target_link_libraries(test_odf_file_io PRIVATE mtrsim Catch2::Catch2WithMain)
target_compile_definitions(test_odf_file_io PRIVATE
  MTRSIM_TEST_DATA_DIR="${MTRSim_SOURCE_DIR}/data"
)
catch_discover_tests(test_odf_file_io)
```

- [ ] **Step 4: Run tests, verify FAIL**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target test_odf_file_io
ctest -R ODFFileIO --output-on-failure
```

**Expected:** Build failure — header exists, no implementation. (If reading tests from an earlier build cache, link error for unresolved `mtrsim::readODFMetadata`.)

- [ ] **Step 5: Implement `readODFMetadata`**

Create `src/LibMTRSim/ODFFileIO.cpp`:

```cpp
#include "ODFFileIO.hpp"

#include <hdf5.h>

#include <stdexcept>
#include <string>

namespace mtrsim {

namespace {

constexpr double k_RadToDeg = 180.0 / 3.14159265358979323846;
constexpr double k_DegToRad = 3.14159265358979323846 / 180.0;

struct FileCloser {
  hid_t id;
  ~FileCloser() { if (id >= 0) H5Fclose(id); }
};

// Read a 1-D dataset of doubles at `path`. Returns size & data.
std::vector<double> readDoubleVector(hid_t fileId, const std::string& path)
{
  const hid_t dset = H5Dopen2(fileId, path.c_str(), H5P_DEFAULT);
  if (dset < 0) {
    throw std::runtime_error("ODFFileIO: cannot open dataset " + path);
  }
  const hid_t space = H5Dget_space(dset);
  hsize_t dims[1] = {0};
  H5Sget_simple_extent_dims(space, dims, nullptr);
  std::vector<double> out(dims[0]);
  if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) < 0) {
    H5Sclose(space); H5Dclose(dset);
    throw std::runtime_error("ODFFileIO: cannot read dataset " + path);
  }
  H5Sclose(space);
  H5Dclose(dset);
  return out;
}

int64_t readScalarInt64(hid_t fileId, const std::string& path)
{
  const hid_t dset = H5Dopen2(fileId, path.c_str(), H5P_DEFAULT);
  if (dset < 0) throw std::runtime_error("ODFFileIO: cannot open scalar " + path);
  int64_t v = 0;
  H5Dread(dset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
  H5Dclose(dset);
  return v;
}

// Validate and derive bin size from an edges array. Returns the uniform step in degrees.
// Throws if edges aren't monotonic/uniform.
double binSizeFromEdgesDeg(const std::vector<double>& edgesRad, const std::string& axisName)
{
  if (edgesRad.size() < 2) {
    throw std::runtime_error("ODFFileIO: " + axisName + "_bins must have >= 2 edges");
  }
  const double step = edgesRad[1] - edgesRad[0];
  constexpr double k_Tol = 1e-12;
  for (std::size_t i = 1; i < edgesRad.size(); ++i) {
    const double d = edgesRad[i] - edgesRad[i - 1];
    if (std::abs(d - step) > k_Tol) {
      throw std::runtime_error("ODFFileIO: " + axisName + "_bins are not uniformly spaced");
    }
    if (d <= 0.0) {
      throw std::runtime_error("ODFFileIO: " + axisName + "_bins are not strictly monotonic");
    }
  }
  return step * k_RadToDeg;
}

}  // namespace

ODFFileMetadata readODFMetadata(const std::filesystem::path& file)
{
  if (!std::filesystem::exists(file)) {
    throw std::runtime_error("ODFFileIO: file does not exist: " + file.string());
  }

  const hid_t fid = H5Fopen(file.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) throw std::runtime_error("ODFFileIO: cannot open " + file.string());
  FileCloser closer{fid};

  ODFFileMetadata meta{};
  meta.numComponents = readScalarInt64(fid, "/ODF_best/num_components");

  if (meta.numComponents < 1) {
    throw std::runtime_error("ODFFileIO: num_components must be >= 1");
  }

  std::vector<double> refPhi1, refPHI, refPhi2;

  for (int64_t c = 0; c < meta.numComponents; ++c) {
    const std::string prefix = "/ODF_best/component_" + std::to_string(c);
    const auto phi1 = readDoubleVector(fid, prefix + "/phi1_bins");
    const auto PHI  = readDoubleVector(fid, prefix + "/PHI_bins");
    const auto phi2 = readDoubleVector(fid, prefix + "/phi2_bins");

    if (c == 0) {
      refPhi1 = phi1; refPHI = PHI; refPhi2 = phi2;
      meta.dimsPhi1PHIPhi2 = {
        static_cast<int64_t>(phi1.size() - 1),
        static_cast<int64_t>(PHI.size()  - 1),
        static_cast<int64_t>(phi2.size() - 1),
      };
      meta.spacingDegPhi1PHIPhi2 = {
        binSizeFromEdgesDeg(phi1, "phi1"),
        binSizeFromEdgesDeg(PHI,  "PHI"),
        binSizeFromEdgesDeg(phi2, "phi2"),
      };
    } else {
      auto bytesMatch = [](const std::vector<double>& a, const std::vector<double>& b) {
        return a.size() == b.size()
          && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
      };
      if (!bytesMatch(phi1, refPhi1) || !bytesMatch(PHI, refPHI) || !bytesMatch(phi2, refPhi2)) {
        throw std::runtime_error(
          "ODFFileIO: component_" + std::to_string(c) + " bin arrays differ from component_0");
      }
    }
  }

  return meta;
}

}  // namespace mtrsim
```

Add `#include <cstring>` for `std::memcmp`.

Also update `src/LibMTRSim/CMakeLists.txt` to include the new source. Find the `target_sources(mtrsim ...)` block and add `ODFFileIO.cpp`.

- [ ] **Step 6: Run tests, verify metadata tests PASS**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target test_odf_file_io
ctest -R ODFFileIO --output-on-failure
```

**Expected:** The two metadata test cases pass. Component-reader and writer tests not yet implemented.

- [ ] **Step 7: Add failing tests for `readODFComponents` and `writeODFFile`**

Append to `tests/test_odf_file_io.cpp`:

```cpp
TEST_CASE("ODFFileIO::readODFComponents returns ODFval payloads", "[ODFFileIO]")
{
  const fs::path odf = testDataDir() / "simulation_ODF.h5";
  const auto components = mtrsim::readODFComponents(odf);
  REQUIRE(components.size() == 3);
  for (const auto& c : components) {
    REQUIRE(c.values.size() == 72 * 36 * 72);  // 186624
  }
  // Sanity: ODF values are non-negative and approximately sum to 1 per component
  for (const auto& c : components) {
    double sum = 0.0;
    for (double v : c.values) {
      REQUIRE(v >= 0.0);
      sum += v;
    }
    REQUIRE(sum == Approx(1.0).epsilon(1e-6));
  }
}

TEST_CASE("ODFFileIO write/read round-trip is byte-exact", "[ODFFileIO]")
{
  const fs::path in  = testDataDir() / "simulation_ODF.h5";
  const fs::path out = fs::temp_directory_path() / "mtrsim_roundtrip.h5";
  if (fs::exists(out)) fs::remove(out);

  const auto meta = mtrsim::readODFMetadata(in);
  const auto comps = mtrsim::readODFComponents(in);
  mtrsim::writeODFFile(out, meta.dimsPhi1PHIPhi2, meta.spacingDegPhi1PHIPhi2, comps);

  const auto meta2  = mtrsim::readODFMetadata(out);
  const auto comps2 = mtrsim::readODFComponents(out);

  REQUIRE(meta2.numComponents == meta.numComponents);
  REQUIRE(meta2.dimsPhi1PHIPhi2 == meta.dimsPhi1PHIPhi2);
  REQUIRE(comps2.size() == comps.size());
  for (std::size_t i = 0; i < comps.size(); ++i) {
    REQUIRE(comps2[i].values == comps[i].values);  // byte-exact Float64 equality
  }
}
```

- [ ] **Step 8: Implement `readODFComponents` and `writeODFFile`**

Append to `src/LibMTRSim/ODFFileIO.cpp`:

```cpp
std::vector<ODFFileComponent> readODFComponents(const std::filesystem::path& file)
{
  // Validate via metadata read first (also checks bin-array consistency).
  const auto meta = readODFMetadata(file);

  const hid_t fid = H5Fopen(file.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) throw std::runtime_error("ODFFileIO: cannot open " + file.string());
  FileCloser closer{fid};

  const std::size_t expectedSize = static_cast<std::size_t>(
    meta.dimsPhi1PHIPhi2[0] * meta.dimsPhi1PHIPhi2[1] * meta.dimsPhi1PHIPhi2[2]);

  std::vector<ODFFileComponent> out;
  out.reserve(meta.numComponents);
  for (int64_t c = 0; c < meta.numComponents; ++c) {
    const std::string path = "/ODF_best/component_" + std::to_string(c) + "/ODFval";
    ODFFileComponent comp;
    comp.values = readDoubleVector(fid, path);
    if (comp.values.size() != expectedSize) {
      throw std::runtime_error(
        "ODFFileIO: component_" + std::to_string(c) + " ODFval size mismatch");
    }
    out.push_back(std::move(comp));
  }
  return out;
}

namespace {

// Generate edge array: edges[i] = i * stepRad for i in [0, N]. Matches MATLAB's
// `0:2*pi/num_bins:2*pi` exactly when stepRad = 2*pi/num_bins.
std::vector<double> makeEdges(int64_t numBins, double stepRad)
{
  std::vector<double> out(static_cast<std::size_t>(numBins + 1));
  for (int64_t i = 0; i <= numBins; ++i) {
    out[static_cast<std::size_t>(i)] = static_cast<double>(i) * stepRad;
  }
  return out;
}

void writeScalarInt64(hid_t fileId, const std::string& path, int64_t v)
{
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t dset = H5Dcreate2(fileId, path.c_str(), H5T_STD_I64LE, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
  H5Dclose(dset);
  H5Sclose(space);
}

void writeDoubleVector(hid_t groupId, const std::string& name, const std::vector<double>& v)
{
  const hsize_t dims[1] = {v.size()};
  const hid_t space = H5Screate_simple(1, dims, nullptr);
  const hid_t dset = H5Dcreate2(groupId, name.c_str(), H5T_IEEE_F64LE, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Dclose(dset);
  H5Sclose(space);
}

}  // namespace

void writeODFFile(const std::filesystem::path& file,
                  const std::array<int64_t, 3>& dims,
                  const std::array<double, 3>& spacingDeg,
                  const std::vector<ODFFileComponent>& components)
{
  if (components.empty()) {
    throw std::runtime_error("ODFFileIO: cannot write file with zero components");
  }
  const std::size_t expectedSize =
    static_cast<std::size_t>(dims[0] * dims[1] * dims[2]);
  for (std::size_t i = 0; i < components.size(); ++i) {
    if (components[i].values.size() != expectedSize) {
      throw std::runtime_error(
        "ODFFileIO: component " + std::to_string(i) + " size does not match dims");
    }
  }

  const hid_t fid = H5Fcreate(file.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) throw std::runtime_error("ODFFileIO: cannot create " + file.string());
  FileCloser closer{fid};

  H5Gclose(H5Gcreate2(fid, "/ODF_best", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeScalarInt64(fid, "/ODF_best/num_components", static_cast<int64_t>(components.size()));

  const auto phi1Edges = makeEdges(dims[0], spacingDeg[0] * k_DegToRad);
  const auto PHIEdges  = makeEdges(dims[1], spacingDeg[1] * k_DegToRad);
  const auto phi2Edges = makeEdges(dims[2], spacingDeg[2] * k_DegToRad);

  for (std::size_t c = 0; c < components.size(); ++c) {
    const std::string groupPath = "/ODF_best/component_" + std::to_string(c);
    const hid_t grp = H5Gcreate2(fid, groupPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    writeDoubleVector(grp, "ODFval",    components[c].values);
    writeDoubleVector(grp, "phi1_bins", phi1Edges);
    writeDoubleVector(grp, "PHI_bins",  PHIEdges);
    writeDoubleVector(grp, "phi2_bins", phi2Edges);
    H5Gclose(grp);
  }
}
```

- [ ] **Step 9: Run all ODFFileIO tests**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target test_odf_file_io
ctest -R ODFFileIO --output-on-failure
```

**Expected:** All four test cases pass. Round-trip test confirms byte-exact preservation.

- [ ] **Step 10: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add src/LibMTRSim/ODFFileIO.{hpp,cpp} src/LibMTRSim/CMakeLists.txt \
        tests/test_odf_file_io.cpp tests/CMakeLists.txt
git commit -m "$(cat <<'EOF'
feat(LibMTRSim): add ODFFileIO helper for MATLAB-format HDF5 round-trip

New mtrsim::readODFMetadata / readODFComponents / writeODFFile API.
Validates byte-exact bin-array consistency across components on read.
Tested against data/simulation_ODF.h5 with byte-exact round-trip.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `ReadMTRSimODFFilter`

**Files:**
- Modify: `src/MTRSim/src/MTRSim/Filters/ReadMTRSimODFFilter.{hpp,cpp}` (from Task 1)
- Modify: `src/MTRSim/src/MTRSim/Filters/Algorithms/ReadMTRSimODF.{hpp,cpp}`
- Modify: `test/ReadMTRSimODFFilterTest.cpp`
- Create: `test/test_data/simulation_ODF_roundtrip.h5.gz` (copy + gzip of `data/simulation_ODF.h5`)

### Context

Reads a MATLAB-format HDF5 file into an `ImageGeom` + N Float64 cell-data arrays, per spec Section 4. Uses `mtrsim::readODFMetadata` in preflight and `mtrsim::readODFComponents` in execute.

SIMPLNX filter pattern: `parameters()` returns a `Parameters` object with typed keys; `preflightImpl()` validates and emits `OutputActions`; `executeImpl()` performs the actual work.

- [ ] **Step 1: Fill in the filter class metadata**

In `ReadMTRSimODFFilter.cpp`, fill in:

```cpp
std::string ReadMTRSimODFFilter::name() const   { return FilterTraits<ReadMTRSimODFFilter>::name; }
std::string ReadMTRSimODFFilter::className() const { return FilterTraits<ReadMTRSimODFFilter>::className; }
Uuid ReadMTRSimODFFilter::uuid() const          { return FilterTraits<ReadMTRSimODFFilter>::uuid; }
std::string ReadMTRSimODFFilter::humanName() const { return "Read MTRSim ODF (HDF5)"; }
std::vector<std::string> ReadMTRSimODFFilter::defaultTags() const {
  return {className(), "IO", "Input", "Read", "Import", "MTRSim", "ODF"};
}
```

- [ ] **Step 2: Define parameters**

Replace the `parameters()` TODO block:

```cpp
Parameters ReadMTRSimODFFilter::parameters() const
{
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<FileSystemPathParameter>(
    k_InputFile_Key, "Input HDF5 File",
    "MATLAB-format MTRSim ODF HDF5 file (.h5 / .hdf5).",
    fs::path{}, FileSystemPathParameter::ExtensionsType{".h5", ".hdf5"},
    FileSystemPathParameter::PathType::InputFile));

  params.insertSeparator(Parameters::Separator{"Output Parameter(s)"});
  params.insert(std::make_unique<DataGroupCreationParameter>(
    k_OutputImageGeometry_Key, "Output ODF Image Geometry",
    "Path to the new ImageGeometry representing the ODF Euler-space grid.",
    DataPath({"ODF"})));
  params.insert(std::make_unique<DataObjectNameParameter>(
    k_CellAttrMatName_Key, "Cell Attribute Matrix Name",
    "Name of the cell attribute matrix that will hold the ODF component arrays.",
    "Cell Data"));

  return params;
}
```

Add matching `inline constexpr StringLiteral` keys near the top of the `.cpp`:

```cpp
namespace {
constexpr StringLiteral k_InputFile_Key          = "input_file";
constexpr StringLiteral k_OutputImageGeometry_Key = "output_image_geometry";
constexpr StringLiteral k_CellAttrMatName_Key    = "cell_attribute_matrix_name";
}
```

Add includes:

```cpp
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "LibMTRSim/ODFFileIO.hpp"
```

- [ ] **Step 3: Implement preflight**

Replace the `preflightImpl` TODO block:

```cpp
IFilter::PreflightResult ReadMTRSimODFFilter::preflightImpl(
  const DataStructure& dataStructure, const Arguments& args,
  const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel) const
{
  const auto inputFile        = args.value<FileSystemPathParameter::ValueType>(k_InputFile_Key);
  const auto outputGeomPath   = args.value<DataGroupCreationParameter::ValueType>(k_OutputImageGeometry_Key);
  const auto cellAttrMatName  = args.value<DataObjectNameParameter::ValueType>(k_CellAttrMatName_Key);

  mtrsim::ODFFileMetadata meta;
  try {
    meta = mtrsim::readODFMetadata(inputFile);
  } catch (const std::exception& e) {
    return {MakeErrorResult<OutputActions>(-12001, std::string("ODF file read failed: ") + e.what())};
  }

  OutputActions actions;

  // ImageGeom axes: (X,Y,Z) = (phi2, PHI, phi1)
  const std::vector<size_t> dimsXYZ = {
    static_cast<size_t>(meta.dimsPhi1PHIPhi2[2]),  // phi2 = X
    static_cast<size_t>(meta.dimsPhi1PHIPhi2[1]),  // PHI  = Y
    static_cast<size_t>(meta.dimsPhi1PHIPhi2[0]),  // phi1 = Z
  };
  const std::vector<float> spacingXYZ = {
    static_cast<float>(meta.spacingDegPhi1PHIPhi2[2]),
    static_cast<float>(meta.spacingDegPhi1PHIPhi2[1]),
    static_cast<float>(meta.spacingDegPhi1PHIPhi2[0]),
  };
  const std::vector<float> originXYZ = {0.0f, 0.0f, 0.0f};

  const DataPath cellAttrMatPath = outputGeomPath.createChildPath(cellAttrMatName);
  actions.appendAction(std::make_unique<CreateImageGeometryAction>(
    outputGeomPath, dimsXYZ, originXYZ, spacingXYZ, cellAttrMatName));

  // Cell data arrays: tuple shape is ZYX-ordered for SIMPLNX. N × single-component Float64.
  const std::vector<size_t> tupleShapeZYX = {dimsXYZ[2], dimsXYZ[1], dimsXYZ[0]};
  const std::vector<size_t> compShape = {1};
  for (int64_t c = 0; c < meta.numComponents; ++c) {
    const std::string arrayName = "component_" + std::to_string(c);
    actions.appendAction(std::make_unique<CreateArrayAction>(
      DataType::float64, tupleShapeZYX, compShape, cellAttrMatPath.createChildPath(arrayName)));
  }

  return {std::move(actions)};
}
```

Add includes:

```cpp
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
```

- [ ] **Step 4: Implement execute**

Replace the `executeImpl` TODO block:

```cpp
Result<> ReadMTRSimODFFilter::executeImpl(
  DataStructure& dataStructure, const Arguments& args,
  const PipelineFilter* /*pipelineNode*/, const MessageHandler& messageHandler,
  const std::atomic_bool& /*shouldCancel*/) const
{
  const auto inputFile       = args.value<FileSystemPathParameter::ValueType>(k_InputFile_Key);
  const auto outputGeomPath  = args.value<DataGroupCreationParameter::ValueType>(k_OutputImageGeometry_Key);
  const auto cellAttrMatName = args.value<DataObjectNameParameter::ValueType>(k_CellAttrMatName_Key);
  const DataPath cellAttrMatPath = outputGeomPath.createChildPath(cellAttrMatName);

  messageHandler({IFilter::Message::Type::Info, "Reading " + inputFile.string()});
  const auto components = mtrsim::readODFComponents(inputFile);

  for (std::size_t c = 0; c < components.size(); ++c) {
    const std::string arrayName = "component_" + std::to_string(c);
    auto& out = dataStructure.getDataRefAs<Float64Array>(
      cellAttrMatPath.createChildPath(arrayName));
    auto& store = out.getDataStoreRef();
    const auto& src = components[c].values;
    if (store.getSize() != src.size()) {
      return MakeErrorResult(-12010,
        "component_" + std::to_string(c) + " tuple count mismatch after preflight");
    }
    for (std::size_t i = 0; i < src.size(); ++i) {
      store.setValue(i, src[i]);
    }
  }
  return {};
}
```

Add include:

```cpp
#include "simplnx/DataStructure/DataArray.hpp"
```

- [ ] **Step 5: Create test exemplar**

```bash
cd /Users/mjackson/Workspace7/MTRSim
mkdir -p test/test_data
cp data/simulation_ODF.h5 test/test_data/
gzip -9 test/test_data/simulation_ODF.h5
ls -lh test/test_data/
```

**Expected:** `simulation_ODF.h5.gz` present (~1 MB).

- [ ] **Step 6: Write the filter-level test**

Replace the scaffolded `ReadMTRSimODFFilterTest.cpp` with:

```cpp
#include "MTRSim/Filters/ReadMTRSimODFFilter.hpp"

#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
using namespace nx::core;

namespace {
// Decompress a .gz fixture from test_data/ into a temp file. Returns the path.
fs::path decompressFixture(const std::string& name)
{
  const fs::path src = fs::path(MTRSIM_TEST_DATA_DIR) / name;
  const fs::path dst = fs::temp_directory_path() / fs::path(name).stem(); // strip .gz
  if (fs::exists(dst)) return dst;
  const std::string cmd = "gunzip -c \"" + src.string() + "\" > \"" + dst.string() + "\"";
  REQUIRE(std::system(cmd.c_str()) == 0);
  return dst;
}
}  // namespace

TEST_CASE("ReadMTRSimODFFilter: round-trip import populates ImageGeom + components",
         "[MTRSim][ReadMTRSimODFFilter]")
{
  const fs::path input = decompressFixture("simulation_ODF.h5.gz");

  ReadMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key,      std::make_any<FileSystemPathParameter::ValueType>(input));
  args.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key, std::make_any<std::string>("Cell Data"));

  DataStructure ds;
  auto pre = filter.preflight(ds, args);
  SIMPLNX_RESULT_REQUIRE_VALID(pre.outputActions);

  auto exec = filter.execute(ds, args);
  SIMPLNX_RESULT_REQUIRE_VALID(exec.result);

  const auto& geom = ds.getDataRefAs<ImageGeom>(DataPath({"ODF"}));
  REQUIRE(geom.getNumXCells() == 72);
  REQUIRE(geom.getNumYCells() == 36);
  REQUIRE(geom.getNumZCells() == 72);
  REQUIRE(geom.getSpacing()[0] == Catch::Approx(5.0));

  for (int c = 0; c < 3; ++c) {
    const std::string name = "component_" + std::to_string(c);
    const auto& arr = ds.getDataRefAs<Float64Array>(
      DataPath({"ODF", "Cell Data", name}));
    REQUIRE(arr.getNumberOfTuples() == 72 * 36 * 72);
    double sum = 0.0;
    for (std::size_t i = 0; i < arr.getNumberOfTuples(); ++i) sum += arr[i];
    REQUIRE(sum == Catch::Approx(1.0).epsilon(1e-6));
  }
}

TEST_CASE("ReadMTRSimODFFilter: missing file errors at preflight",
         "[MTRSim][ReadMTRSimODFFilter][ErrorPath]")
{
  ReadMTRSimODFFilter filter;
  Arguments args;
  args.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key,
                     std::make_any<FileSystemPathParameter::ValueType>(fs::path("/nonexistent.h5")));
  args.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key,
                     std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key,
                     std::make_any<std::string>("Cell Data"));

  DataStructure ds;
  auto pre = filter.preflight(ds, args);
  REQUIRE(pre.outputActions.invalid());
}
```

Update the test's CMakeLists (auto-generated by `make_filter.py`) to pass `MTRSIM_TEST_DATA_DIR`:

```cmake
target_compile_definitions(MTRSimTest PRIVATE
  MTRSIM_TEST_DATA_DIR="${${PLUGIN_NAME}_SOURCE_DIR}/test/test_data"
)
```

- [ ] **Step 7: Run filter tests**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target MTRSimTest
ctest -R ReadMTRSimODF --output-on-failure
```

**Expected:** Both test cases pass.

- [ ] **Step 8: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add src/MTRSim/src/MTRSim/Filters/ReadMTRSimODFFilter.{hpp,cpp} \
        src/MTRSim/src/MTRSim/Filters/Algorithms/ReadMTRSimODF.{hpp,cpp} \
        test/ReadMTRSimODFFilterTest.cpp \
        test/test_data/simulation_ODF.h5.gz
git commit -m "$(cat <<'EOF'
feat(MTRSim): implement ReadMTRSimODFFilter

Reads a MATLAB-format MTRSim ODF HDF5 file into an ImageGeom with N
Float64 cell-data arrays (one per component). Preflight uses
mtrsim::readODFMetadata; execute uses readODFComponents. Tested
against a gzipped exemplar.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `WriteMTRSimODFFilter` (+ round-trip test)

**Files:**
- Modify: `src/MTRSim/src/MTRSim/Filters/WriteMTRSimODFFilter.{hpp,cpp}`
- Modify: `src/MTRSim/src/MTRSim/Filters/Algorithms/WriteMTRSimODF.{hpp,cpp}`
- Modify: `test/WriteMTRSimODFFilterTest.cpp`

### Context

Writes a SIMPLNX ODF (ImageGeom + selected Float64 cell arrays) back out to the MATLAB-compatible HDF5 format, per spec Section 5. Round-trips losslessly with `ReadMTRSimODFFilter`. Uses `mtrsim::writeODFFile` from Task 2.

- [ ] **Step 1: Fill in class metadata + parameters**

In `WriteMTRSimODFFilter.cpp`, fill in `name()` / `humanName()` / `defaultTags()` analogously to Task 3 (human name: `"Write MTRSim ODF (HDF5)"`; tags include `"Output", "Write", "Export"`). Add parameters:

```cpp
Parameters WriteMTRSimODFFilter::parameters() const
{
  Parameters params;
  params.insertSeparator(Parameters::Separator{"Output Parameter(s)"});
  params.insert(std::make_unique<FileSystemPathParameter>(
    k_OutputFile_Key, "Output HDF5 File",
    "Destination file for the MATLAB-format MTRSim ODF.",
    fs::path{}, FileSystemPathParameter::ExtensionsType{".h5", ".hdf5"},
    FileSystemPathParameter::PathType::OutputFile));

  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<GeometrySelectionParameter>(
    k_InputImageGeometry_Key, "Input ODF Image Geometry",
    "ImageGeom representing the Euler-space ODF grid.",
    DataPath{},
    GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<MultiArraySelectionParameter>(
    k_ODFComponents_Key, "ODF Components to Export",
    "Cell-data arrays on the selected geometry to serialize as component_0..component_{N-1}.",
    MultiArraySelectionParameter::ValueType{},
    MultiArraySelectionParameter::AllowedTypes{IArray::ArrayType::DataArray},
    MultiArraySelectionParameter::AllowedDataTypes{DataType::float64},
    MultiArraySelectionParameter::AllowedComponentShapes{{1}}));

  return params;
}
```

With key constants:

```cpp
namespace {
constexpr StringLiteral k_OutputFile_Key          = "output_file";
constexpr StringLiteral k_InputImageGeometry_Key  = "input_image_geometry";
constexpr StringLiteral k_ODFComponents_Key       = "odf_components";
}
```

- [ ] **Step 2: Implement preflight**

```cpp
IFilter::PreflightResult WriteMTRSimODFFilter::preflightImpl(
  const DataStructure& dataStructure, const Arguments& args, ...) const
{
  const auto geomPath = args.value<DataPath>(k_InputImageGeometry_Key);
  const auto arrayPaths = args.value<MultiArraySelectionParameter::ValueType>(k_ODFComponents_Key);

  if (arrayPaths.empty()) {
    return {MakeErrorResult<OutputActions>(-12100, "Select at least one ODF component array.")};
  }

  const auto* geom = dataStructure.getDataAs<ImageGeom>(geomPath);
  if (geom == nullptr) {
    return {MakeErrorResult<OutputActions>(-12101, "Input geometry must be an ImageGeom.")};
  }

  const size_t expectedTuples =
    geom->getNumXCells() * geom->getNumYCells() * geom->getNumZCells();

  for (const auto& p : arrayPaths) {
    // Ensure the array is a descendant of the geometry's cell attribute matrix
    if (!p.getParent().toString().starts_with(geomPath.toString())) {
      return {MakeErrorResult<OutputActions>(-12102,
        "Array '" + p.toString() + "' is not on the selected geometry.")};
    }
    const auto* arr = dataStructure.getDataAs<Float64Array>(p);
    if (arr == nullptr) {
      return {MakeErrorResult<OutputActions>(-12103,
        "Array '" + p.toString() + "' must be Float64.")};
    }
    if (arr->getNumberOfTuples() != expectedTuples) {
      return {MakeErrorResult<OutputActions>(-12104,
        "Array tuple count does not match geometry cell count.")};
    }
  }
  if (geom->getNumXCells() < 2 || geom->getNumYCells() < 2 || geom->getNumZCells() < 2) {
    return {MakeErrorResult<OutputActions>(-12105, "Degenerate geometry rejected (< 2 bins on some axis).")};
  }
  return {OutputActions{}};
}
```

- [ ] **Step 3: Implement execute**

```cpp
Result<> WriteMTRSimODFFilter::executeImpl(
  DataStructure& dataStructure, const Arguments& args, ...) const
{
  const auto outputFile  = args.value<FileSystemPathParameter::ValueType>(k_OutputFile_Key);
  const auto geomPath    = args.value<DataPath>(k_InputImageGeometry_Key);
  const auto arrayPaths  = args.value<MultiArraySelectionParameter::ValueType>(k_ODFComponents_Key);

  const auto& geom = dataStructure.getDataRefAs<ImageGeom>(geomPath);
  const auto spacingXYZ = geom.getSpacing();

  // Remap geometry (X=phi2, Y=PHI, Z=phi1) back to (phi1, PHI, phi2) for the file format.
  const std::array<int64_t, 3> dimsPhi1PHIPhi2 = {
    static_cast<int64_t>(geom.getNumZCells()),
    static_cast<int64_t>(geom.getNumYCells()),
    static_cast<int64_t>(geom.getNumXCells()),
  };
  const std::array<double, 3> spacingDegPhi1PHIPhi2 = {
    static_cast<double>(spacingXYZ[2]),
    static_cast<double>(spacingXYZ[1]),
    static_cast<double>(spacingXYZ[0]),
  };

  std::vector<mtrsim::ODFFileComponent> components;
  components.reserve(arrayPaths.size());
  for (const auto& p : arrayPaths) {
    const auto& arr = dataStructure.getDataRefAs<Float64Array>(p);
    mtrsim::ODFFileComponent comp;
    comp.values.resize(arr.getNumberOfTuples());
    for (std::size_t i = 0; i < comp.values.size(); ++i) comp.values[i] = arr[i];
    components.push_back(std::move(comp));
  }

  try {
    mtrsim::writeODFFile(outputFile, dimsPhi1PHIPhi2, spacingDegPhi1PHIPhi2, components);
  } catch (const std::exception& e) {
    return MakeErrorResult(-12110, std::string("ODF write failed: ") + e.what());
  }
  return {};
}
```

- [ ] **Step 4: Write round-trip test**

In `WriteMTRSimODFFilterTest.cpp`:

```cpp
TEST_CASE("WriteMTRSimODFFilter: Import→Export produces byte-exact round-trip",
         "[MTRSim][WriteMTRSimODFFilter]")
{
  const fs::path input  = decompressFixture("simulation_ODF.h5.gz");
  const fs::path output = fs::temp_directory_path() / "mtrsim_export_roundtrip.h5";
  if (fs::exists(output)) fs::remove(output);

  DataStructure ds;
  // --- Import
  {
    ReadMTRSimODFFilter imp;
    Arguments a;
    a.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(input));
    a.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));
    a.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key, std::make_any<std::string>("Cell Data"));
    SIMPLNX_RESULT_REQUIRE_VALID(imp.preflight(ds, a).outputActions);
    SIMPLNX_RESULT_REQUIRE_VALID(imp.execute(ds, a).result);
  }
  // --- Export
  {
    WriteMTRSimODFFilter exp;
    Arguments a;
    a.insertOrAssign(WriteMTRSimODFFilter::k_OutputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(output));
    a.insertOrAssign(WriteMTRSimODFFilter::k_InputImageGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));
    std::vector<DataPath> comps = {
      DataPath({"ODF","Cell Data","component_0"}),
      DataPath({"ODF","Cell Data","component_1"}),
      DataPath({"ODF","Cell Data","component_2"}),
    };
    a.insertOrAssign(WriteMTRSimODFFilter::k_ODFComponents_Key, std::make_any<std::vector<DataPath>>(comps));
    SIMPLNX_RESULT_REQUIRE_VALID(exp.preflight(ds, a).outputActions);
    SIMPLNX_RESULT_REQUIRE_VALID(exp.execute(ds, a).result);
  }

  // --- Verify: byte-exact on ODFval arrays (bin edges will match because step is identical)
  const auto origComps    = mtrsim::readODFComponents(input);
  const auto roundTripped = mtrsim::readODFComponents(output);
  REQUIRE(origComps.size() == roundTripped.size());
  for (std::size_t i = 0; i < origComps.size(); ++i) {
    REQUIRE(origComps[i].values == roundTripped[i].values);
  }
}

TEST_CASE("WriteMTRSimODFFilter: rejects zero component selection",
         "[MTRSim][WriteMTRSimODFFilter][ErrorPath]")
{
  DataStructure ds;
  // Build a minimal ImageGeom by hand (or skip — just test preflight)
  WriteMTRSimODFFilter exp;
  Arguments a;
  a.insertOrAssign(WriteMTRSimODFFilter::k_OutputFile_Key,
                  std::make_any<FileSystemPathParameter::ValueType>(fs::temp_directory_path() / "unused.h5"));
  a.insertOrAssign(WriteMTRSimODFFilter::k_InputImageGeometry_Key, std::make_any<DataPath>(DataPath{}));
  a.insertOrAssign(WriteMTRSimODFFilter::k_ODFComponents_Key, std::make_any<std::vector<DataPath>>({}));
  auto pre = exp.preflight(ds, a);
  REQUIRE(pre.outputActions.invalid());
}
```

- [ ] **Step 5: Run tests**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target MTRSimTest
ctest -R WriteMTRSimODF --output-on-failure
```

**Expected:** Round-trip + error-path tests pass.

- [ ] **Step 6: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add src/MTRSim/src/MTRSim/Filters/WriteMTRSimODFFilter.{hpp,cpp} \
        src/MTRSim/src/MTRSim/Filters/Algorithms/WriteMTRSimODF.{hpp,cpp} \
        test/WriteMTRSimODFFilterTest.cpp
git commit -m "$(cat <<'EOF'
feat(MTRSim): implement WriteMTRSimODFFilter

Writes selected Float64 cell-data arrays on an ODF ImageGeom to the
MATLAB-format HDF5 layout. Round-trip test confirms byte-exact
preservation of ODFval arrays through Import -> Export -> read-back.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: `mtrsim::SymmetricEulers` library helper (TDD)

**Files:**
- Create: `src/LibMTRSim/SymmetricEulers.{hpp,cpp}`
- Create: `tests/test_symmetric_eulers.cpp`
- Modify: `src/LibMTRSim/CMakeLists.txt`
- Modify: `tests/CMakeLists.txt`

### Context

Wraps EbsdLib's `LaueOps::getMatSymOpD(i)` API to produce all symmetrically-equivalent Euler tuples for a given Bunge orientation + crystal-system code. Analog of MATLAB's `symmetric_euler_angles.m`, but routed through EbsdLib so orientation math stays centralized. Pure library code; no simplnx dependency. Independent of Task 2 — good for parallel subagent execution.

- [ ] **Step 1: Write the header**

Create `src/LibMTRSim/SymmetricEulers.hpp`:

```cpp
#pragma once

#include "libmtrsim_export.h"

#include <array>
#include <cstdint>
#include <vector>

namespace mtrsim {

/// Returns the list of symmetrically-equivalent Bunge Euler tuples (phi1, PHI, phi2)
/// in radians. The input is a single Bunge tuple and a crystal-structure code
/// matching EbsdLib's EbsdLib::CrystalStructure:: enum (e.g., Hexagonal_High = 0).
/// The returned list length depends on the crystal system (e.g., 24 for Hexagonal_High).
MTRSIM_EXPORT std::vector<std::array<double, 3>>
expandSymmetric(double phi1Rad, double PHIRad, double phi2Rad, uint32_t ebsdLibCrystalCode);

}  // namespace mtrsim
```

- [ ] **Step 2: Write failing tests**

Create `tests/test_symmetric_eulers.cpp`:

```cpp
#include "LibMTRSim/SymmetricEulers.hpp"

#include "EbsdLib/Core/EbsdLibConstants.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>

using Catch::Approx;

TEST_CASE("SymmetricEulers: HCP input returns 24 equivalents (12 symmetry ops x 2 for inversion)",
         "[SymmetricEulers]")
{
  const auto out = mtrsim::expandSymmetric(0.1, 0.2, 0.3, EbsdLib::CrystalStructure::Hexagonal_High);
  // HCP Laue group (6/mmm) has 24 rotations including improper via inversion handling in LaueOps
  REQUIRE(out.size() == 24);
  // The identity (or its equivalent) must reproduce the input within tolerance on at least one entry
  bool foundIdentity = false;
  for (const auto& e : out) {
    if (std::abs(e[0] - 0.1) < 1e-9
     && std::abs(e[1] - 0.2) < 1e-9
     && std::abs(e[2] - 0.3) < 1e-9) {
      foundIdentity = true; break;
    }
  }
  REQUIRE(foundIdentity);
}

TEST_CASE("SymmetricEulers: cubic returns 24 equivalents", "[SymmetricEulers]")
{
  const auto out = mtrsim::expandSymmetric(0.1, 0.2, 0.3, EbsdLib::CrystalStructure::Cubic_High);
  REQUIRE(out.size() == 24);
}

TEST_CASE("SymmetricEulers: all outputs are valid Bunge angles (in the 0..2pi, 0..pi, 0..2pi fundamental zone)",
         "[SymmetricEulers]")
{
  const auto out = mtrsim::expandSymmetric(0.1, 0.2, 0.3, EbsdLib::CrystalStructure::Hexagonal_High);
  constexpr double k_TwoPi = 2.0 * 3.14159265358979323846;
  constexpr double k_Pi    = 3.14159265358979323846;
  for (const auto& e : out) {
    REQUIRE(e[0] >= 0.0);
    REQUIRE(e[0] <  k_TwoPi + 1e-9);
    REQUIRE(e[1] >= 0.0);
    REQUIRE(e[1] <= k_Pi    + 1e-9);
    REQUIRE(e[2] >= 0.0);
    REQUIRE(e[2] <  k_TwoPi + 1e-9);
  }
}
```

Register in `tests/CMakeLists.txt`:

```cmake
add_executable(test_symmetric_eulers test_symmetric_eulers.cpp)
target_link_libraries(test_symmetric_eulers PRIVATE mtrsim EbsdLib::EbsdLib Catch2::Catch2WithMain)
catch_discover_tests(test_symmetric_eulers)
```

- [ ] **Step 3: Run tests, verify FAIL**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target test_symmetric_eulers 2>&1 | tail -10
```

**Expected:** Link failure — `expandSymmetric` undefined.

- [ ] **Step 4: Implement using EbsdLib**

Create `src/LibMTRSim/SymmetricEulers.cpp`:

```cpp
#include "SymmetricEulers.hpp"

#include "EbsdLib/LaueOps/LaueOps.h"
#include "EbsdLib/Core/EbsdLibConstants.h"
#include "EbsdLib/Core/OrientationTransformation.hpp"

#include <stdexcept>

namespace mtrsim {

std::vector<std::array<double, 3>>
expandSymmetric(double phi1Rad, double PHIRad, double phi2Rad, uint32_t ebsdLibCrystalCode)
{
  using OT = EbsdLib::OrientationTransformation;

  // Convert Bunge -> rotation matrix
  const std::array<double, 3> eu = {phi1Rad, PHIRad, phi2Rad};
  const std::array<double, 9> om = OT::eu2om<std::array<double, 3>, std::array<double, 9>>(eu);

  // Get LaueOps subclass for crystal system
  auto ops = LaueOps::GetLaueOpsFromCrystalStructure(ebsdLibCrystalCode);
  if (ops == nullptr) {
    throw std::runtime_error("SymmetricEulers: unknown crystal structure code "
                           + std::to_string(ebsdLibCrystalCode));
  }
  const size_t nOps = ops->getNumSymOps();

  std::vector<std::array<double, 3>> out;
  out.reserve(nOps);

  for (size_t i = 0; i < nOps; ++i) {
    const auto sym = ops->getMatSymOpD(i);  // Matrix3X3D

    // Multiply: om_out = sym * om
    std::array<double, 9> omOut{};
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c)
        for (int k = 0; k < 3; ++k)
          omOut[r*3 + c] += sym(r, k) * om[k*3 + c];

    const auto euOut = OT::om2eu<std::array<double, 9>, std::array<double, 3>>(omOut);
    out.push_back({euOut[0], euOut[1], euOut[2]});
  }
  return out;
}

}  // namespace mtrsim
```

Add to `src/LibMTRSim/CMakeLists.txt`:

```cmake
target_sources(mtrsim PRIVATE SymmetricEulers.cpp)
```

(Ensure `EbsdLib::EbsdLib` is a PUBLIC link target of `mtrsim` — per MEMORY.md it already is.)

- [ ] **Step 5: Run tests, verify PASS**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target test_symmetric_eulers
ctest -R SymmetricEulers --output-on-failure
```

**Expected:** All three test cases pass. If the identity-check fails, likely because `getMatSymOpD(0)` isn't the identity on all subclasses — relax the test to look for an orientation that reproduces `(phi1, PHI, phi2)` within tolerance (angle equivalence allows for 2π periodicity in phi1/phi2).

- [ ] **Step 6: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add src/LibMTRSim/SymmetricEulers.{hpp,cpp} src/LibMTRSim/CMakeLists.txt \
        tests/test_symmetric_eulers.cpp tests/CMakeLists.txt
git commit -m "$(cat <<'EOF'
feat(LibMTRSim): add SymmetricEulers helper (EbsdLib-backed)

Returns all symmetrically-equivalent Bunge Euler tuples for an input
orientation + crystal structure. Wraps EbsdLib::LaueOps symmetry
operators, keeping all orientation math centralized in EbsdLib.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: `mtrsim::ODFBuilder` library helper (TDD)

**Files:**
- Create: `src/LibMTRSim/ODFBuilder.{hpp,cpp}`
- Create: `tests/test_odf_builder.cpp`
- Modify: `src/LibMTRSim/CMakeLists.txt`
- Modify: `tests/CMakeLists.txt`

### Context

The core binning + smoothing accumulator. Given a list of symmetric Euler tuples (from `SymmetricEulers::expandSymmetric`) and an `ODFBuildParams`, accumulate contributions into an ODFval vector with tri-linear face/edge/corner weights matching `matlab/calc_ODF.m`. Depends on Task 5.

- [ ] **Step 1: Write the header**

Create `src/LibMTRSim/ODFBuilder.hpp`:

```cpp
#pragma once

#include "libmtrsim_export.h"

#include <array>
#include <cstdint>
#include <vector>

namespace mtrsim {

struct MTRSIM_EXPORT ODFBuildParams {
  int32_t nphi1;         // bins along phi1 (slowest-varying)
  int32_t nPHI;          // bins along PHI
  int32_t nphi2;         // bins along phi2 (fastest-varying)
  double  binSizeDeg;    // uniform across all three axes
  bool    smoothing;     // tri-linear face/edge/corner distribution when true
};

/// Accumulates contributions from a single set of symmetric Euler tuples into
/// `values`. Safe to call repeatedly; does NOT normalize. `values` must be
/// pre-sized to nphi1 * nPHI * nphi2 and initially zero.
/// Eulers are expected in radians.
MTRSIM_EXPORT void accumulate(const std::vector<std::array<double, 3>>& symmetricEulers,
                              const ODFBuildParams& params,
                              std::vector<double>& values);

/// Final normalization: values /= normalizer.
MTRSIM_EXPORT void normalize(std::vector<double>& values, double normalizer);

}  // namespace mtrsim
```

- [ ] **Step 2: Write failing tests**

Create `tests/test_odf_builder.cpp`:

```cpp
#include "LibMTRSim/ODFBuilder.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <numbers>

using Catch::Approx;

TEST_CASE("ODFBuilder::accumulate (no smoothing) hits exactly one bin per tuple", "[ODFBuilder]")
{
  constexpr double k_DegToRad = std::numbers::pi / 180.0;
  const mtrsim::ODFBuildParams p{
    .nphi1 = 72, .nPHI = 36, .nphi2 = 72, .binSizeDeg = 5.0, .smoothing = false
  };
  std::vector<double> values(72 * 36 * 72, 0.0);

  // A single tuple at bin (i_phi1=5, i_PHI=3, i_phi2=10) center: angles (27.5, 17.5, 52.5) deg
  const std::vector<std::array<double, 3>> eulers = {
    {27.5 * k_DegToRad, 17.5 * k_DegToRad, 52.5 * k_DegToRad}
  };
  mtrsim::accumulate(eulers, p, values);

  // Expect 1.0 at the target bin, 0.0 elsewhere
  double total = 0.0;
  for (double v : values) total += v;
  REQUIRE(total == Approx(1.0).epsilon(1e-12));

  const std::size_t ix = static_cast<std::size_t>(5 * 36 * 72 + 3 * 72 + 10);
  REQUIRE(values[ix] == Approx(1.0).epsilon(1e-12));
}

TEST_CASE("ODFBuilder::accumulate (smoothing) distributes per MATLAB weights", "[ODFBuilder]")
{
  constexpr double k_DegToRad = std::numbers::pi / 180.0;
  const mtrsim::ODFBuildParams p{
    .nphi1 = 72, .nPHI = 36, .nphi2 = 72, .binSizeDeg = 5.0, .smoothing = true
  };
  std::vector<double> values(72 * 36 * 72, 0.0);

  const std::vector<std::array<double, 3>> eulers = {
    {27.5 * k_DegToRad, 17.5 * k_DegToRad, 52.5 * k_DegToRad}
  };
  mtrsim::accumulate(eulers, p, values);

  // Total contribution per tuple = 0.332 + 0.448 + 0.16 + 0.06 = 1.0
  double total = 0.0;
  for (double v : values) total += v;
  REQUIRE(total == Approx(1.0).epsilon(1e-9));

  // Center bin should hold exactly 0.332
  const std::size_t ix = static_cast<std::size_t>(5 * 36 * 72 + 3 * 72 + 10);
  REQUIRE(values[ix] == Approx(0.332).epsilon(1e-12));
}

TEST_CASE("ODFBuilder::normalize divides in place", "[ODFBuilder]")
{
  std::vector<double> v = {2.0, 4.0, 8.0};
  mtrsim::normalize(v, 2.0);
  REQUIRE(v == std::vector<double>{1.0, 2.0, 4.0});
}
```

Register in `tests/CMakeLists.txt`:

```cmake
add_executable(test_odf_builder test_odf_builder.cpp)
target_link_libraries(test_odf_builder PRIVATE mtrsim Catch2::Catch2WithMain)
catch_discover_tests(test_odf_builder)
```

- [ ] **Step 3: Run tests, verify FAIL**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target test_odf_builder 2>&1 | tail -5
```

**Expected:** Link error.

- [ ] **Step 4: Implement**

Create `src/LibMTRSim/ODFBuilder.cpp`:

```cpp
#include "ODFBuilder.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace mtrsim {

namespace {

constexpr double k_CenterWeight  = 0.332;
constexpr double k_FaceWeight    = 0.448 / 6.0;
constexpr double k_EdgeWeight    = 0.16  / 12.0;
constexpr double k_CornerWeight  = 0.06  / 8.0;

inline int32_t wrap(int32_t i, int32_t n) { return ((i % n) + n) % n; }

inline std::size_t linearize(int32_t i_phi1, int32_t i_PHI, int32_t i_phi2,
                             int32_t nphi1, int32_t nPHI, int32_t nphi2)
{
  return static_cast<std::size_t>(
    i_phi1 * (nPHI * nphi2) + i_PHI * nphi2 + i_phi2);
}

}  // namespace

void accumulate(const std::vector<std::array<double, 3>>& symmetricEulers,
                const ODFBuildParams& p,
                std::vector<double>& values)
{
  const std::size_t total =
    static_cast<std::size_t>(p.nphi1) * p.nPHI * p.nphi2;
  if (values.size() != total) {
    throw std::invalid_argument("ODFBuilder::accumulate: values size mismatch");
  }
  const double k_RadToDeg = 180.0 / std::numbers::pi;

  for (const auto& eu : symmetricEulers) {
    const double phi1Deg = eu[0] * k_RadToDeg;
    const double PHIDeg  = eu[1] * k_RadToDeg;
    const double phi2Deg = eu[2] * k_RadToDeg;

    const int32_t i_phi1 = wrap(static_cast<int32_t>(std::floor(phi1Deg / p.binSizeDeg)), p.nphi1);
    const int32_t i_PHI  = wrap(static_cast<int32_t>(std::floor(PHIDeg  / p.binSizeDeg)), p.nPHI);
    const int32_t i_phi2 = wrap(static_cast<int32_t>(std::floor(phi2Deg / p.binSizeDeg)), p.nphi2);

    if (!p.smoothing) {
      values[linearize(i_phi1, i_PHI, i_phi2, p.nphi1, p.nPHI, p.nphi2)] += 1.0;
      continue;
    }

    // Center
    values[linearize(i_phi1, i_PHI, i_phi2, p.nphi1, p.nPHI, p.nphi2)] += k_CenterWeight;

    // Six face neighbors (offset ±1 along exactly one axis)
    constexpr int faceOffsets[6][3] = {
      {-1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, +1, 0}, {0, 0, -1}, {0, 0, +1}
    };
    for (auto& o : faceOffsets) {
      values[linearize(
        wrap(i_phi1 + o[0], p.nphi1),
        wrap(i_PHI  + o[1], p.nPHI),
        wrap(i_phi2 + o[2], p.nphi2),
        p.nphi1, p.nPHI, p.nphi2)] += k_FaceWeight;
    }

    // Twelve edge neighbors (offset ±1 along two axes)
    constexpr int edgeOffsets[12][3] = {
      {-1,-1, 0}, {-1,+1, 0}, {+1,-1, 0}, {+1,+1, 0},
      {-1, 0,-1}, {-1, 0,+1}, {+1, 0,-1}, {+1, 0,+1},
      { 0,-1,-1}, { 0,-1,+1}, { 0,+1,-1}, { 0,+1,+1}
    };
    for (auto& o : edgeOffsets) {
      values[linearize(
        wrap(i_phi1 + o[0], p.nphi1),
        wrap(i_PHI  + o[1], p.nPHI),
        wrap(i_phi2 + o[2], p.nphi2),
        p.nphi1, p.nPHI, p.nphi2)] += k_EdgeWeight;
    }

    // Eight corner neighbors (offset ±1 along all three axes)
    for (int sp1 : {-1, +1}) for (int sP : {-1, +1}) for (int sp2 : {-1, +1}) {
      values[linearize(
        wrap(i_phi1 + sp1, p.nphi1),
        wrap(i_PHI  + sP,  p.nPHI),
        wrap(i_phi2 + sp2, p.nphi2),
        p.nphi1, p.nPHI, p.nphi2)] += k_CornerWeight;
    }
  }
}

void normalize(std::vector<double>& values, double normalizer)
{
  if (normalizer == 0.0) return;
  for (double& v : values) v /= normalizer;
}

}  // namespace mtrsim
```

Add to `src/LibMTRSim/CMakeLists.txt`: `target_sources(mtrsim PRIVATE ODFBuilder.cpp)`.

- [ ] **Step 5: Run tests, verify PASS**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target test_odf_builder
ctest -R ODFBuilder --output-on-failure
```

**Expected:** All three test cases pass.

- [ ] **Step 6: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add src/LibMTRSim/ODFBuilder.{hpp,cpp} src/LibMTRSim/CMakeLists.txt \
        tests/test_odf_builder.cpp tests/CMakeLists.txt
git commit -m "$(cat <<'EOF'
feat(LibMTRSim): add ODFBuilder with MATLAB tri-linear smoothing weights

Accumulates Euler-tuple contributions into a row-major ODFval vector,
with optional tri-linear face/edge/corner smoothing using the
MATLAB calc_ODF.m weights (0.332 / 0.448/6 / 0.16/12 / 0.06/8). Wrap
boundaries use Bunge-angle periodicity modulo the bin count.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: `ComputeODFFilter`

**Files:**
- Modify: `src/MTRSim/src/MTRSim/Filters/ComputeODFFilter.{hpp,cpp}`
- Modify: `src/MTRSim/src/MTRSim/Filters/Algorithms/ComputeODF.{hpp,cpp}`
- Modify: `test/ComputeODFFilterTest.cpp`
- Create: `test/test_data/hcp_euler_sample.h5.gz` (EBSD fixture, produced below)
- Create: `test/test_data/hcp_euler_sample_odf_reference.h5.gz` (MATLAB-generated reference)

### Context

Spec Section 6. Consumes Euler angles + phases + crystal structures (+ optional mask), uses `mtrsim::expandSymmetric` (Task 5) and `mtrsim::ODFBuilder::accumulate / normalize` (Task 6) to build the ODF. Supports both "Create New" and "Append to Existing" modes.

- [ ] **Step 1: Fill in class metadata + parameters**

In `ComputeODFFilter.cpp`:

```cpp
std::string ComputeODFFilter::humanName() const { return "Compute ODF From Euler Angles"; }
std::vector<std::string> ComputeODFFilter::defaultTags() const {
  return {className(), "ODF", "MTRSim", "Orientation", "Statistics"};
}

Parameters ComputeODFFilter::parameters() const
{
  Parameters params;

  params.insertLinkableParameter(std::make_unique<ChoicesParameter>(
    k_OutputMode_Key, "Output Mode",
    "Create a new ODF Image Geometry or append a new component to an existing one.",
    0, ChoicesParameter::Choices{"Create New ODF Geometry", "Append to Existing ODF Geometry"}));

  params.insertSeparator(Parameters::Separator{"Algorithm Parameters"});
  params.insert(std::make_unique<BoolParameter>(
    k_Smoothing_Key, "Apply Tri-linear Smoothing", "If true, distribute each tuple's contribution "
    "across the center bin and 26 neighbors per MATLAB calc_ODF.m.", true));
  params.insert(std::make_unique<Float32Parameter>(
    k_BinSizeDeg_Key, "Bin Size (degrees)",
    "Uniform across phi1/PHI/phi2. Must divide 360 (phi1, phi2) and 180 (PHI) evenly.", 5.0f));

  params.insertSeparator(Parameters::Separator{"Input Arrays"});
  params.insert(std::make_unique<ArraySelectionParameter>(
    k_EulerAngles_Key, "Euler Angles",
    "Float32, 3-component, per-cell Bunge Euler angles (radians).",
    DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::float32},
    ArraySelectionParameter::AllowedComponentShapes{{3}}));
  params.insert(std::make_unique<ArraySelectionParameter>(
    k_Phases_Key, "Phases",
    "Int32, single-component, per-cell phase indices.",
    DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::int32},
    ArraySelectionParameter::AllowedComponentShapes{{1}}));
  params.insert(std::make_unique<ArraySelectionParameter>(
    k_CrystalStructures_Key, "Crystal Structures",
    "UInt32, single-component, per-phase (ensemble-level) crystal structure codes.",
    DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::uint32},
    ArraySelectionParameter::AllowedComponentShapes{{1}}));
  params.insertLinkableParameter(std::make_unique<BoolParameter>(
    k_UseMask_Key, "Use Mask Array", "Limit contributions to masked-true voxels.", false));
  params.insert(std::make_unique<ArraySelectionParameter>(
    k_Mask_Key, "Mask Array",
    "Bool, single-component, per-cell mask.",
    DataPath{}, ArraySelectionParameter::AllowedTypes{DataType::boolean},
    ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insertSeparator(Parameters::Separator{"Create New Mode"});
  params.insert(std::make_unique<DataGroupCreationParameter>(
    k_NewGeom_Key, "New ODF Image Geometry", "Path for the newly created ImageGeom.", DataPath({"ODF"})));
  params.insert(std::make_unique<DataObjectNameParameter>(
    k_NewCellAttrMat_Key, "Cell Attribute Matrix Name", "", "Cell Data"));

  params.insertSeparator(Parameters::Separator{"Append Mode"});
  params.insert(std::make_unique<GeometrySelectionParameter>(
    k_ExistingGeom_Key, "Existing ODF Image Geometry",
    "Target for the new component array.",
    DataPath{},
    GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<DataObjectNameParameter>(
    k_ExistingCellAttrMat_Key, "Existing Cell Attribute Matrix Name",
    "Name of the existing cell attribute matrix on the selected geometry where the new array will live.",
    "Cell Data"));

  params.insert(std::make_unique<DataObjectNameParameter>(
    k_ComponentName_Key, "New Component Array Name",
    "Name of the new Float64 cell-data array.", "Component 1"));

  params.linkParameters(k_OutputMode_Key, k_NewGeom_Key,             std::make_any<ChoicesParameter::ValueType>(0));
  params.linkParameters(k_OutputMode_Key, k_NewCellAttrMat_Key,      std::make_any<ChoicesParameter::ValueType>(0));
  params.linkParameters(k_OutputMode_Key, k_BinSizeDeg_Key,          std::make_any<ChoicesParameter::ValueType>(0));
  params.linkParameters(k_OutputMode_Key, k_ExistingGeom_Key,        std::make_any<ChoicesParameter::ValueType>(1));
  params.linkParameters(k_OutputMode_Key, k_ExistingCellAttrMat_Key, std::make_any<ChoicesParameter::ValueType>(1));
  params.linkParameters(k_UseMask_Key,    k_Mask_Key,                std::make_any<bool>(true));
  return params;
}
```

Key constants:

```cpp
namespace {
constexpr StringLiteral k_OutputMode_Key        = "output_mode";
constexpr StringLiteral k_Smoothing_Key         = "apply_smoothing";
constexpr StringLiteral k_BinSizeDeg_Key        = "bin_size_deg";
constexpr StringLiteral k_EulerAngles_Key       = "euler_angles";
constexpr StringLiteral k_Phases_Key            = "phases";
constexpr StringLiteral k_CrystalStructures_Key = "crystal_structures";
constexpr StringLiteral k_UseMask_Key           = "use_mask";
constexpr StringLiteral k_Mask_Key              = "mask";
constexpr StringLiteral k_NewGeom_Key           = "new_geom";
constexpr StringLiteral k_NewCellAttrMat_Key    = "new_cell_attr_mat";
constexpr StringLiteral k_ExistingGeom_Key        = "existing_geom";
constexpr StringLiteral k_ExistingCellAttrMat_Key = "existing_cell_attr_mat";
constexpr StringLiteral k_ComponentName_Key       = "component_name";
}
```

- [ ] **Step 2: Implement preflight**

```cpp
IFilter::PreflightResult ComputeODFFilter::preflightImpl(
  const DataStructure& ds, const Arguments& args, ...) const
{
  const auto mode       = args.value<ChoicesParameter::ValueType>(k_OutputMode_Key);
  const float binDeg    = args.value<float>(k_BinSizeDeg_Key);
  const auto eulerPath  = args.value<DataPath>(k_EulerAngles_Key);
  const auto phasesPath = args.value<DataPath>(k_Phases_Key);
  const auto componentName = args.value<std::string>(k_ComponentName_Key);

  const auto* eulers = ds.getDataAs<Float32Array>(eulerPath);
  const auto* phases = ds.getDataAs<Int32Array>(phasesPath);
  if (eulers == nullptr || phases == nullptr) {
    return {MakeErrorResult<OutputActions>(-12200, "Euler angles and phases arrays must exist.")};
  }
  if (phases->getNumberOfTuples() != eulers->getNumberOfTuples()) {
    return {MakeErrorResult<OutputActions>(-12201, "Euler and phases tuple counts differ.")};
  }

  OutputActions actions;
  size_t nphi1 = 0, nPHI = 0, nphi2 = 0;
  DataPath arrayPath;

  if (mode == 0) {  // Create New
    if (binDeg <= 0.0f) return {MakeErrorResult<OutputActions>(-12202, "bin_size_deg must be > 0")};
    if (std::fmod(360.0f, binDeg) > 1e-4f || std::fmod(180.0f, binDeg) > 1e-4f) {
      return {MakeErrorResult<OutputActions>(-12203,
        "bin_size_deg must divide 360 (phi1, phi2) and 180 (PHI) evenly.")};
    }
    nphi2 = static_cast<size_t>(std::round(360.0f / binDeg));
    nPHI  = static_cast<size_t>(std::round(180.0f / binDeg));
    nphi1 = static_cast<size_t>(std::round(360.0f / binDeg));

    const auto newGeom = args.value<DataPath>(k_NewGeom_Key);
    const auto cellMat = args.value<std::string>(k_NewCellAttrMat_Key);
    actions.appendAction(std::make_unique<CreateImageGeometryAction>(
      newGeom, std::vector<size_t>{nphi2, nPHI, nphi1},
      std::vector<float>{0.0f, 0.0f, 0.0f},
      std::vector<float>{binDeg, binDeg, binDeg}, cellMat));
    arrayPath = newGeom.createChildPath(cellMat).createChildPath(componentName);
  } else {  // Append
    const auto existing = args.value<DataPath>(k_ExistingGeom_Key);
    const auto existingCellMat = args.value<std::string>(k_ExistingCellAttrMat_Key);
    const auto* geom = ds.getDataAs<ImageGeom>(existing);
    if (geom == nullptr) return {MakeErrorResult<OutputActions>(-12210, "Existing geometry must be ImageGeom.")};
    const auto sp = geom->getSpacing();
    if (std::abs(sp[0] - sp[1]) > 1e-6f || std::abs(sp[1] - sp[2]) > 1e-6f) {
      return {MakeErrorResult<OutputActions>(-12211, "Existing geometry spacing must be uniform across axes.")};
    }
    nphi2 = geom->getNumXCells();
    nPHI  = geom->getNumYCells();
    nphi1 = geom->getNumZCells();
    // User explicitly names the cell-attr-matrix to avoid brittle ImageGeom API lookups.
    const DataPath cellMatPath = existing.createChildPath(existingCellMat);
    if (ds.getDataAs<AttributeMatrix>(cellMatPath) == nullptr) {
      return {MakeErrorResult<OutputActions>(-12213,
        "Cell attribute matrix '" + cellMatPath.toString() + "' not found on the existing geometry.")};
    }
    arrayPath = cellMatPath.createChildPath(componentName);

    // Collision check
    if (ds.getDataAs<IDataArray>(arrayPath) != nullptr) {
      return {MakeErrorResult<OutputActions>(-12212,
        "Array '" + arrayPath.toString() + "' already exists on the target geometry.")};
    }
  }

  const std::vector<size_t> tupleShapeZYX = {nphi1, nPHI, nphi2};
  actions.appendAction(std::make_unique<CreateArrayAction>(
    DataType::float64, tupleShapeZYX, std::vector<size_t>{1}, arrayPath));

  return {std::move(actions)};
}
```

- [ ] **Step 3: Implement execute**

Pull most of the logic into the Algorithm class (`ComputeODF`):

```cpp
// ComputeODF.hpp
namespace nx::core::MTRSim {

struct ComputeODFInputValues {
  DataPath eulerAnglesPath;
  DataPath phasesPath;
  DataPath crystalStructuresPath;
  DataPath maskPath;            // empty if not used
  DataPath outputArrayPath;
  int32_t  nphi1, nPHI, nphi2;
  double   binSizeDeg;
  bool     smoothing;
};

class ComputeODF
{
public:
  ComputeODF(DataStructure&, const ComputeODFInputValues*,
                           const IFilter::MessageHandler&, const std::atomic_bool&);
  Result<> operator()();
private:
  // ...
};

}  // namespace
```

Algorithm body (in `.cpp`):

```cpp
Result<> ComputeODF::operator()()
{
  const auto& eulers = m_DataStructure.getDataRefAs<Float32Array>(m_Inputs->eulerAnglesPath);
  const auto& phases = m_DataStructure.getDataRefAs<Int32Array>(m_Inputs->phasesPath);
  const auto& xtal   = m_DataStructure.getDataRefAs<UInt32Array>(m_Inputs->crystalStructuresPath);
  const BoolArray* mask = m_Inputs->maskPath.empty()
    ? nullptr
    : m_DataStructure.getDataAs<BoolArray>(m_Inputs->maskPath);
  auto& out = m_DataStructure.getDataRefAs<Float64Array>(m_Inputs->outputArrayPath);

  const mtrsim::ODFBuildParams p{
    m_Inputs->nphi1, m_Inputs->nPHI, m_Inputs->nphi2,
    m_Inputs->binSizeDeg, m_Inputs->smoothing
  };

  std::vector<double> values(
    static_cast<std::size_t>(p.nphi1) * p.nPHI * p.nphi2, 0.0);

  const std::size_t N = eulers.getNumberOfTuples();
  std::size_t contributing = 0;
  for (std::size_t i = 0; i < N; ++i) {
    if (mask != nullptr && !(*mask)[i]) continue;
    const int32_t phaseId = phases[i];
    if (phaseId < 0 || static_cast<std::size_t>(phaseId) >= xtal.getNumberOfTuples()) continue;
    if (phaseId == 0) continue;  // unindexed
    const uint32_t xs = xtal[static_cast<std::size_t>(phaseId)];
    const double phi1 = eulers[i * 3 + 0];
    const double PHI  = eulers[i * 3 + 1];
    const double phi2 = eulers[i * 3 + 2];
    const auto sym = mtrsim::expandSymmetric(phi1, PHI, phi2, xs);
    mtrsim::accumulate(sym, p, values);
    ++contributing;
  }
  if (contributing == 0) return MakeErrorResult(-12250, "No voxels contributed to the ODF.");
  mtrsim::normalize(values, static_cast<double>(contributing));

  for (std::size_t i = 0; i < values.size(); ++i) out[i] = values[i];
  return {};
}
```

The filter's `executeImpl` just populates the input-values struct and invokes the algorithm class.

- [ ] **Step 4: Prepare test fixtures**

Create a small HCP fixture and its MATLAB reference. This is done **once**, committed, and reused:

```bash
cd /Users/mjackson/Workspace7/MTRSim
# 1. Extract a small HCP sample from existing data if available, or build synthetic one
#    via a throw-away helper script. 100-500 Euler tuples is sufficient for test purposes.
# 2. Run matlab/calc_ODF.m on this sample to produce the reference ODFval.
# 3. Save both to HDF5 and gzip:
#      hcp_euler_sample.h5             — {phi1, PHI, phi2 arrays; num tuples; crystal code}
#      hcp_euler_sample_odf_reference.h5 — {ODFval array, binSizeDeg, smoothing flag}
gzip -9 test/test_data/hcp_euler_sample.h5
gzip -9 test/test_data/hcp_euler_sample_odf_reference.h5
```

**Note:** If MATLAB is not immediately available on the implementation machine, defer this step until the test-data phase, stub the test with `WARN("reference missing; skipping tolerance check")`, and file a follow-up task. The algorithmic tests in Task 6 already verify the binning/smoothing unit logic against known numeric weights.

- [ ] **Step 5: Write the filter-level test**

Minimum viable test in `ComputeODFFilterTest.cpp`:

```cpp
// Fixture schema (hcp_euler_sample.h5):
//   /num_points            int64                 — N
//   /crystal_code          uint32                — EbsdLib::CrystalStructure code
//   /phi1                  float64[N]            — radians
//   /PHI                   float64[N]            — radians
//   /phi2                  float64[N]            — radians
// Fixture schema (hcp_euler_sample_odf_reference.h5):
//   /bin_size_deg          float64 scalar
//   /smoothing             uint8   scalar (0 or 1)
//   /odf_val               float64[nphi1*nPHI*nphi2]  — MATLAB calc_ODF.m reference

TEST_CASE("ComputeODFFilter: Create New mode against MATLAB reference",
         "[MTRSim][ComputeODFFilter]")
{
  const fs::path euler = decompressFixture("hcp_euler_sample.h5.gz");
  const fs::path ref   = decompressFixture("hcp_euler_sample_odf_reference.h5.gz");

  // Read fixture
  hid_t fid = H5Fopen(euler.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  int64_t N = 0; uint32_t crystalCode = 0;
  // [use H5Dread to fill N, crystalCode, phi1/PHI/phi2 vectors]
  int64_t Nlocal = N; uint32_t codeLocal = crystalCode;
  std::vector<double> phi1Vec(Nlocal), PHIVec(Nlocal), phi2Vec(Nlocal);
  // ... read datasets ...
  H5Fclose(fid);

  // Read reference
  hid_t rid = H5Fopen(ref.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  double binSizeDeg = 5.0; uint8_t smoothing = 1;
  std::vector<double> refODFVal(72 * 36 * 72);
  // ... read datasets ...
  H5Fclose(rid);

  // Build DataStructure with the required input arrays on a tiny source ImageGeom.
  DataStructure ds;
  {
    auto* srcGeom = ImageGeom::Create(ds, "EBSD");
    srcGeom->setDimensions({static_cast<size_t>(Nlocal), 1, 1});
    srcGeom->setSpacing({1.0f, 1.0f, 1.0f});
    srcGeom->setOrigin({0.0f, 0.0f, 0.0f});
    // ensemble matrix with crystal_structures
    auto* ensMat = AttributeMatrix::Create(ds, "Ensemble", {2}, srcGeom->getId());
    auto* xtals = Float64Array::CreateWithStore<DataStore<uint32_t>>(ds, "CrystalStructures", {2}, {1}, ensMat->getId());
    (*xtals)[0] = 999;               // unindexed sentinel
    (*xtals)[1] = codeLocal;
    // cell matrix with eulers, phases
    auto* cellMat = AttributeMatrix::Create(ds, "Cell Data", {static_cast<size_t>(Nlocal)}, srcGeom->getId());
    auto* eulers = Float32Array::CreateWithStore<DataStore<float>>(ds, "EulerAngles",
                                                                 {static_cast<size_t>(Nlocal)}, {3}, cellMat->getId());
    auto* phases = Int32Array::CreateWithStore<DataStore<int32_t>>(ds, "Phases",
                                                                 {static_cast<size_t>(Nlocal)}, {1}, cellMat->getId());
    for (int64_t i = 0; i < Nlocal; ++i) {
      (*eulers)[i*3 + 0] = static_cast<float>(phi1Vec[i]);
      (*eulers)[i*3 + 1] = static_cast<float>(PHIVec[i]);
      (*eulers)[i*3 + 2] = static_cast<float>(phi2Vec[i]);
      (*phases)[i] = 1;
    }
  }

  ComputeODFFilter filter;
  Arguments args;
  args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(0));
  args.insertOrAssign(ComputeODFFilter::k_Smoothing_Key,  std::make_any<bool>(smoothing != 0));
  args.insertOrAssign(ComputeODFFilter::k_BinSizeDeg_Key, std::make_any<float>(static_cast<float>(binSizeDeg)));
  args.insertOrAssign(ComputeODFFilter::k_EulerAngles_Key,
                     std::make_any<DataPath>(DataPath({"EBSD","Cell Data","EulerAngles"})));
  args.insertOrAssign(ComputeODFFilter::k_Phases_Key,
                     std::make_any<DataPath>(DataPath({"EBSD","Cell Data","Phases"})));
  args.insertOrAssign(ComputeODFFilter::k_CrystalStructures_Key,
                     std::make_any<DataPath>(DataPath({"EBSD","Ensemble","CrystalStructures"})));
  args.insertOrAssign(ComputeODFFilter::k_UseMask_Key, std::make_any<bool>(false));
  args.insertOrAssign(ComputeODFFilter::k_NewGeom_Key, std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(ComputeODFFilter::k_NewCellAttrMat_Key, std::make_any<std::string>("Cell Data"));
  args.insertOrAssign(ComputeODFFilter::k_ComponentName_Key, std::make_any<std::string>("Component 1"));

  SIMPLNX_RESULT_REQUIRE_VALID(filter.preflight(ds, args).outputActions);
  SIMPLNX_RESULT_REQUIRE_VALID(filter.execute(ds, args).result);

  const auto& out = ds.getDataRefAs<Float64Array>(DataPath({"ODF","Cell Data","Component 1"}));
  REQUIRE(out.getNumberOfTuples() == refODFVal.size());
  double maxAbsDiff = 0.0;
  for (std::size_t i = 0; i < refODFVal.size(); ++i) {
    maxAbsDiff = std::max(maxAbsDiff, std::abs(out[i] - refODFVal[i]));
  }
  INFO("maxAbsDiff = " << maxAbsDiff);
  REQUIRE(maxAbsDiff < 1e-10);
}

TEST_CASE("ComputeODFFilter: Append mode adds a new component",
         "[MTRSim][ComputeODFFilter]")
{
  // Seed with an ODF geometry produced by ReadMTRSimODFFilter, then run
  // ComputeODF in Append mode with the same fixture. After execute, the
  // target geometry should carry N+1 components where N is the original count.
  const fs::path input = decompressFixture("simulation_ODF.h5.gz");
  const fs::path euler = decompressFixture("hcp_euler_sample.h5.gz");

  DataStructure ds;
  // Import first to create /ODF geometry with component_0..component_2
  {
    ReadMTRSimODFFilter imp;
    Arguments a;
    a.insertOrAssign(ReadMTRSimODFFilter::k_InputFile_Key,           std::make_any<FileSystemPathParameter::ValueType>(input));
    a.insertOrAssign(ReadMTRSimODFFilter::k_OutputImageGeometry_Key, std::make_any<DataPath>(DataPath({"ODF"})));
    a.insertOrAssign(ReadMTRSimODFFilter::k_CellAttrMatName_Key,     std::make_any<std::string>("Cell Data"));
    SIMPLNX_RESULT_REQUIRE_VALID(imp.preflight(ds, a).outputActions);
    SIMPLNX_RESULT_REQUIRE_VALID(imp.execute(ds, a).result);
  }
  // [Populate EBSD arrays per the fixture, same as Create-New test above]

  ComputeODFFilter filter;
  Arguments args;
  args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(1));
  args.insertOrAssign(ComputeODFFilter::k_ExistingGeom_Key,
                     std::make_any<DataPath>(DataPath({"ODF"})));
  args.insertOrAssign(ComputeODFFilter::k_ExistingCellAttrMat_Key,
                     std::make_any<std::string>("Cell Data"));
  args.insertOrAssign(ComputeODFFilter::k_ComponentName_Key,
                     std::make_any<std::string>("Component 3"));
  // [plus the common input/smoothing args — per Create-New test]

  SIMPLNX_RESULT_REQUIRE_VALID(filter.preflight(ds, args).outputActions);
  SIMPLNX_RESULT_REQUIRE_VALID(filter.execute(ds, args).result);

  const auto& newArr = ds.getDataRefAs<Float64Array>(DataPath({"ODF","Cell Data","Component 3"}));
  REQUIRE(newArr.getNumberOfTuples() == 72 * 36 * 72);
  // Original three components remain untouched
  REQUIRE(ds.getDataAs<Float64Array>(DataPath({"ODF","Cell Data","component_0"})) != nullptr);
  REQUIRE(ds.getDataAs<Float64Array>(DataPath({"ODF","Cell Data","component_2"})) != nullptr);
}

TEST_CASE("ComputeODFFilter: non-integer bin size rejected at preflight",
         "[MTRSim][ComputeODFFilter][ErrorPath]")
{
  DataStructure ds;
  // Build the minimum ensemble+cell arrays required to pass the earlier input checks.
  // [Populate dummy ensemble + single-tuple cell-level arrays as in Create-New test]

  ComputeODFFilter filter;
  Arguments args;
  args.insertOrAssign(ComputeODFFilter::k_OutputMode_Key, std::make_any<ChoicesParameter::ValueType>(0));
  args.insertOrAssign(ComputeODFFilter::k_BinSizeDeg_Key, std::make_any<float>(7.0f));  // 180/7 non-integer
  // [plus the common input args]

  auto pre = filter.preflight(ds, args);
  REQUIRE(pre.outputActions.invalid());
}
```

- [ ] **Step 6: Build + run tests**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake --build . --target MTRSimTest
ctest -R ComputeODF --output-on-failure
```

**Expected:** Happy path passes within tolerance; error-path rejects non-integer bin size.

- [ ] **Step 7: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add src/MTRSim/src/MTRSim/Filters/ComputeODFFilter.{hpp,cpp} \
        src/MTRSim/src/MTRSim/Filters/Algorithms/ComputeODF.{hpp,cpp} \
        test/ComputeODFFilterTest.cpp \
        test/test_data/hcp_euler_sample.h5.gz \
        test/test_data/hcp_euler_sample_odf_reference.h5.gz
git commit -m "$(cat <<'EOF'
feat(MTRSim): implement ComputeODFFilter

Builds an ODF from per-voxel Bunge Euler angles + phase data, with
optional tri-linear smoothing. Uses EbsdLib-backed symmetry expansion
via mtrsim::expandSymmetric, and mtrsim::ODFBuilder for binning.
Supports Create New and Append to Existing output modes. Tolerance-
based test against a MATLAB calc_ODF.m reference fixture.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Filter documentation

**Files:**
- Modify: `docs/ReadMTRSimODFFilter.md`
- Modify: `docs/WriteMTRSimODFFilter.md`
- Modify: `docs/ComputeODFFilter.md`

### Context

Follow the BlueQuartz `filter-documentation` skill's template. Each filter gets one markdown file covering purpose, parameters, data requirements, output, and an example use case. The axis-mapping note from spec Section 3 must appear in each filter's doc.

- [ ] **Step 1: Write ReadMTRSimODFFilter.md**

Follow the template:

```markdown
# Read MTRSim ODF (HDF5)

## Group (Subgroup)

MTRSim (Input)

## Description

Reads a MATLAB-format MTRSim ODF HDF5 file into an ImageGeometry + N Float64 cell-data arrays.

## Axis Mapping

**Important:** phi1 is the slowest-varying axis and maps to the ImageGeometry **Z** dimension.
PHI maps to Y, phi2 maps to X. To view standard ODF plots (phi1–PHI plane at fixed phi2),
use DREAM3D-NX's YZ image-slicing view at a chosen X slice.

## Parameters

| Name | Type | Description |
|------|------|-------------|
| Input HDF5 File | File | `.h5` or `.hdf5` containing `/ODF_best/num_components` + `component_N` groups |
| Output ODF Image Geometry | DataPath | Where the new ImageGeom is created (default `ODF`) |
| Cell Attribute Matrix Name | string | Name of the cell attribute matrix (default `Cell Data`) |

## Required Objects

None — the input file is the sole source.

## Created Objects

- **Image Geometry** at the selected path. Spacing in degrees (phi1/PHI/phi2 bin size).
- **N cell-data arrays** named `component_0 .. component_{N-1}`, one per ODF component, Float64.

## Example Pipeline

1. Read MTRSim ODF (HDF5) → creates `/ODF`.
2. Compute IPF Colors (downstream, for visualization).

## Authors

BlueQuartz Software, LLC

## License

Commercial license (DREAM3D-NX plugin).
```

- [ ] **Step 2: Write WriteMTRSimODFFilter.md**

```markdown
# Write MTRSim ODF (HDF5)

## Group (Subgroup)

MTRSim (Output)

## Description

Writes an ODF ImageGeometry's selected Float64 cell-data arrays to a MATLAB-compatible HDF5
file. Produces the same layout documented in `app/main.cpp:60-65`. Preserves byte-exact
Float64 values; regenerates `phi1_bins`, `PHI_bins`, `phi2_bins` edge arrays from the
geometry's spacing.

## Axis Mapping

Internally the ImageGeom maps (X, Y, Z) to (phi2, PHI, phi1). On export this remapping is
reversed so the file layout matches the MATLAB convention (phi1 slowest-varying).

## Parameters

| Name | Type | Description |
|------|------|-------------|
| Output HDF5 File | File | Destination `.h5` / `.hdf5` path |
| Input ODF Image Geometry | DataPath | ImageGeom holding the ODF components |
| ODF Components to Export | MultiArraySelection | Float64, single-component cell-data arrays to serialize |

## Semantics

Exported HDF5 groups are named `component_0..component_{N-1}` where `i` is the **selection
index in this filter's array list**, not the source array's data path. Users wanting
control over which array maps to which HDF5 index must select in the desired order. This
is documented and not enforced.

## Required Objects

- Image Geometry with at least one selected Float64 cell-data array.

## Created Objects

- One HDF5 file on disk.

## Example Pipeline

1. Import MTRSim ODF → `/ODF`.
2. (Optionally modify) compute something, rename arrays, etc.
3. Export MTRSim ODF — select all `component_*` arrays; write `round_trip.h5`.

## Authors

BlueQuartz Software, LLC
```

- [ ] **Step 3: Write ComputeODFFilter.md**

```markdown
# Compute ODF From Euler Angles

## Group (Subgroup)

MTRSim (Orientation)

## Description

Builds an Orientation Distribution Function (ODF) from per-voxel Bunge Euler angles and
phase data. Applies crystal-symmetry expansion (via EbsdLib's LaueOps) and optional
tri-linear neighbor-bin smoothing matching the weights in MATLAB's `calc_ODF.m`. Output
is a Float64 cell-data array on an ImageGeom Euler-space grid — compatible with
`WriteMTRSimODFFilter`.

## Axis Mapping

ImageGeom: (X, Y, Z) = (phi2, PHI, phi1). The 5° default uses dimensions 72 × 36 × 72 =
186,624 bins.

## Modes

- **Create New ODF Geometry**: creates a new ImageGeom + cell attribute matrix + the
  first component array.
- **Append to Existing ODF Geometry**: adds a new named Float64 component to an existing
  ODF ImageGeom. Used for multi-phase or multi-scan accumulation workflows (run the
  filter once per phase/scan, each producing a named component).

## Smoothing Weights

When **Apply Tri-linear Smoothing** is true, each symmetric Euler tuple deposits 1.0 of
total contribution distributed as follows (matching `matlab/calc_ODF.m`):

| Target | Count | Per-bin weight |
|--------|-------|----------------|
| Center bin | 1 | 0.332 |
| Face neighbors | 6 | 0.448 / 6 |
| Edge neighbors | 12 | 0.16 / 12 |
| Corner neighbors | 8 | 0.06 / 8 |

Boundaries wrap under Bunge-angle periodicity.

## Parameters

| Name | Type | Description |
|------|------|-------------|
| Output Mode | Choices | Create New ODF Geometry \| Append to Existing ODF Geometry |
| Apply Tri-linear Smoothing | Bool | Default: on |
| Bin Size (degrees) | Float32 | Uniform. Must divide 360 (phi1, phi2) and 180 (PHI) evenly. Hidden/ignored in Append mode |
| Euler Angles | ArraySelection | Float32, 3-component, cell-level, Bunge radians |
| Phases | ArraySelection | Int32, single-component, cell-level |
| Crystal Structures | ArraySelection | UInt32, single-component, ensemble-level. Values match `EbsdLib::CrystalStructure` codes |
| Use Mask Array | Bool | Toggles mask visibility |
| Mask Array | ArraySelection | Bool, cell-level. Only masked-true voxels contribute |
| (Create New) New ODF Image Geometry | DataGroupCreation | Default: `ODF` |
| (Create New) Cell Attribute Matrix Name | string | Default: `Cell Data` |
| (Append) Existing ODF Image Geometry | GeometrySelection | ImageGeom to receive the new component |
| (Append) Existing Cell Attribute Matrix Name | string | Default: `Cell Data` |
| New Component Array Name | string | Default: `Component 1` |

## Required Objects

- Cell-level Float32 3-component Euler-angles array (radians).
- Cell-level Int32 phases array.
- Ensemble-level UInt32 crystal-structures array (usually produced by the **Create
  Ensemble Info** filter).
- (Optional) Cell-level Bool mask array.

## Created Objects

- (Create New) An ImageGeom + cell attribute matrix + one Float64 single-component array.
- (Append) One Float64 single-component array on the selected existing geometry.

## Lineage

Port of MATLAB `calc_ODF.m` + `symmetric_euler_angles.m`. Symmetry operators now route
through EbsdLib (`LaueOps::getMatSymOpD`) rather than a hand-rolled HCP matrix list.

## Example Pipeline

1. Read H5Ebsd → `/DataContainer`.
2. Create Ensemble Info → adds `CrystalStructures` to the ensemble matrix.
3. Threshold Array (optional) → build a mask restricting to α-Ti voxels.
4. **Compute ODF From Euler Angles** (Create New mode) → creates `/ODF`.
5. Export MTRSim ODF → write a MATLAB-compatible file for validation.

## Authors

BlueQuartz Software, LLC
```

- [ ] **Step 4: Verify `DOC_CHECK` passes during configure**

```bash
cd /Users/mjackson/Workspace7/Build/mtrsim-Rel
cmake . 2>&1 | grep -i "doc"
```

**Expected:** No errors about missing filter documentation.

- [ ] **Step 5: Commit**

```bash
cd /Users/mjackson/Workspace7/MTRSim
git add docs/
git commit -m "$(cat <<'EOF'
docs(MTRSim): filter documentation for Import/Export/ComputeODF

Adds user-facing documentation for the three AJ filters, covering
parameters, required/created objects, axis-mapping caveat, and
example-pipeline snippets.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Deferred decisions to revisit

- **Exemplar storage location** (Task 3 currently stores gzipped HDF5 in-repo under `test/test_data/`). Revisit when a fixture exceeds ~10 MB; candidates are GitHub release artifacts and a shared simplnx webserver (see spec Section 10).
- **Tolerance threshold** for the ComputeODF filter's MATLAB-reference test (Task 7 uses 1e-10 tentatively). Tune once the first comparison is run.
- **Append-mode `bin_size_deg` visibility** — Task 7 currently uses `linkParameters` to hide it; if the simplnx convention for grayed-out-but-visible differs, switch accordingly.

---

## Self-review checklist (run before executing)

- [ ] Every spec requirement maps to at least one task — check Section 4 (Import), 5 (Export), 6 (Compute), 7 (Testing), 8 (Acceptance), 9 (Out-of-scope) against tasks above.
- [ ] No placeholder text — only actual code or concrete commands.
- [ ] Type consistency — `mtrsim::ODFFileMetadata`, `ODFBuildParams.nphi1/nPHI/nphi2`, `expandSymmetric` signature used consistently across tasks.
- [ ] Filter human-names, tags, and UUIDs from `make_filter.py` (UUID is auto-generated; don't hand-write).

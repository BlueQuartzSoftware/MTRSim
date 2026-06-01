# Milestone AK — `MTRSimFilter` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver the `MTRSimFilter` ("Generate Synthetic Microtexture") DREAM3D-NX filter that runs the full MTR plurigaussian simulation from a pipeline, writing MTR ids, Euler angles, and optional polar-color cell data into a new Image Geometry, backed by deterministic + statistical unit tests and green CI.

**Architecture:** The heavy lifting is consolidated into a single reusable LibMTRSim entry point `simulateMTR()` that runs PGRF → orientation sampling → per-voxel assignment and returns results in **SIMPLNX z,y,x voxel order**. LibMTRSim Catch2 tests (in `tests/`) cover the numerics (deterministic helpers + a statistical end-to-end test reusing the ODF exemplar). The filter is a thin SIMPLNX wrapper (`MTRSimFilter` + `MTRSim` algorithm) whose tests (in `test/`) only verify the filter's value-add: parameter→params mapping, ODF-geometry→component reconstruction, array creation/types/names/1-based ids, seed handling, preflight validation, and the optional color array.

**Tech Stack:** C++17, Eigen, simplnx filter framework, Catch2 v2 (note: `Catch::Approx`, not `Approx`), CMake, GitHub Actions.

---

## Background the engineer needs

**Voxel orderings (the highest-risk detail):**
- `GPGenerator::generate` and therefore `PGRFResult::mtrIndex` are ordered
  `kSim = iz·(nx·ny) + ix·ny + iy` (z slowest, **y fastest** — MATLAB
  column-major for an `ny×nx×nz` array). See `src/LibMTRSim/GPGenerator.cpp:53-117`.
- SIMPLNX Image Geometry cell arrays must be ordered
  `kNx = iz·(ny·nx) + iy·nx + ix` (z slowest, **x fastest**).
- The remap from sim order to SIMPLNX order is therefore:
  `out[ iz·(ny·nx) + iy·nx + ix ] = in[ iz·(nx·ny) + ix·ny + iy ]`.

**Grid dimensions:** `nx = round(xLen/dx)`, `ny = round(yLen/dy)`,
`nz = max(round(zLen/dz), 1)`. For the filter: `xLen=Size[0]`, `yLen=Size[1]`,
`zLen=Size[2]`; `dx=Spacing[0]` etc.

**Existing patterns to mirror exactly:**
- Filter boilerplate: `src/MTRSim/Filters/ReadMTRSimODFFilter.{hpp,cpp}`.
- Algorithm boilerplate: `src/MTRSim/Filters/Algorithms/ReadMTRSimODF.{hpp,cpp}`.
- Random-seed parameter pattern: `simplnx` `MergeTwinsFilter.cpp` (lines 60-63,
  103, 184-191).
- ImageGeom + per-array preflight actions: `ReadMTRSimODFFilter::preflightImpl`.

**Reference ODF reconstruction:** the AJ read path writes ODFval row-major in
ZYX (`ReadMTRSimODFFilter.cpp:124-128`). The filter's ODF→component helper is the
inverse: flat values + grid dims + degree-spacing → `mtrsim::ODFComponent` with
bin centres in radians.

**Build/test commands** (dual build per the simplnx convention; the local
simplnx checkout is at `/Users/mjackson/Workspace7/simplnx`):
- LibMTRSim Catch2 tests build via the standalone preset and run with `ctest`.
  Standalone configure/build: `cmake --preset <local-standalone-preset>` then
  `cmake --build <build>`; run `ctest --test-dir <build> -R MTRSim` (see
  `CMakePresets.json` / `CMakeUserPresets.json`).
- Plugin + filter tests build inside the simplnx tree:
  `cmake --build /Users/mjackson/Workspace7/simplnx/build` then
  `ctest --test-dir /Users/mjackson/Workspace7/simplnx/build -R MTRSim --output-on-failure`.

> Before each "run the test" step, prefer the standalone build for Task 1-4
> (fast, no SIMPLNX) and the simplnx build for Task 5-11.

---

## File Structure

**Create:**
- `src/LibMTRSim/MTRSimDriver.hpp` — `MTRSimResult`, `simulateMTR()`,
  `buildUniformODF()`, `gridToODFComponent()`, `remapSimToZYX()` declarations.
- `src/LibMTRSim/MTRSimDriver.cpp` — implementations (logic lifted from
  `src/app/main.cpp`).
- `src/MTRSim/Filters/MTRSimFilter.{hpp,cpp}` — the filter.
- `src/MTRSim/Filters/Algorithms/MTRSim.{hpp,cpp}` — the algorithm.
- `tests/test_mtrsim_driver.cpp` — Catch2 tests for the driver + helpers.
- `test/MTRSimTest.cpp` — SIMPLNX filter tests.
- `docs/MTRSim/MTRSimFilter.md` — filter documentation stub.

**Modify:**
- `src/app/main.cpp` — call `simulateMTR()` instead of inline orchestration.
- `MTRSimPlugin.cmake` — add `MTRSimDriver.{hpp,cpp}` to LibMTRSim sources;
  add `MTRSimFilter`/`MTRSim` to `FilterList`/`AlgorithmList`.
- `LibMTRSim.cmake` — add `MTRSimDriver.{hpp,cpp}` to the standalone library.
- `tests/CMakeLists.txt` — add `test_mtrsim_driver.cpp`.
- `test/CMakeLists.txt` — add `MTRSimTest.cpp`.

> **No git worktree** for this repo (DREAM3DNX build expects MTRSim at a fixed
> path). Work directly on branch `topic/create_mtr_sim_filter`.

---

## Phase A — LibMTRSim building blocks (fast, deterministic, no SIMPLNX)

### Task 1: `buildUniformODF(n1, nPHI, n2)` in the driver header

**Files:**
- Create: `src/LibMTRSim/MTRSimDriver.hpp`
- Create: `src/LibMTRSim/MTRSimDriver.cpp`
- Modify: `LibMTRSim.cmake`, `MTRSimPlugin.cmake`
- Create: `tests/test_mtrsim_driver.cpp`
- Modify: `tests/CMakeLists.txt`

- [ ] **Step 1: Create the driver header skeleton with `buildUniformODF`**

Create `src/LibMTRSim/MTRSimDriver.hpp`:

```cpp
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

} // namespace mtrsim
```

Create `src/LibMTRSim/MTRSimDriver.cpp` (move the body from
`src/app/main.cpp:118-148`, parameterized by dims):

```cpp
#include "MTRSimDriver.hpp"

#include <numbers>

namespace mtrsim {

ODFComponent buildUniformODF(int n1, int nPHI, int n2) {
  const int nTotal = n1 * nPHI * n2;
  const double twoPiOverN1 = 2.0 * std::numbers::pi / static_cast<double>(n1);
  const double piOverNPHI = std::numbers::pi / static_cast<double>(nPHI);
  const double twoPiOverN2 = 2.0 * std::numbers::pi / static_cast<double>(n2);

  Eigen::VectorXd phi1Bins(nTotal);
  Eigen::VectorXd phiBins(nTotal);
  Eigen::VectorXd phi2Bins(nTotal);

  for (int ix = 0; ix < nTotal; ++ix) {
    const int i1 = ix / (nPHI * n2);
    const int iPHI = (ix % (nPHI * n2)) / n2;
    const int i2 = ix % n2;
    phi1Bins[ix] = (i1 + 0.5) * twoPiOverN1;
    phiBins[ix] = (iPHI + 0.5) * piOverNPHI;
    phi2Bins[ix] = (i2 + 0.5) * twoPiOverN2;
  }

  ODFComponent uni;
  uni.odfVal = Eigen::VectorXd::Constant(nTotal, 1.0 / static_cast<double>(nTotal));
  uni.phi1Bins = std::move(phi1Bins);
  uni.phiBins = std::move(phiBins);
  uni.phi2Bins = std::move(phi2Bins);
  return uni;
}

} // namespace mtrsim
```

- [ ] **Step 2: Register the new files in both CMake lists**

In `LibMTRSim.cmake`, add `MTRSimDriver.cpp`/`.hpp` to the standalone library's
source/header lists (mirror the existing `ODFSampler.cpp` entries).

In `MTRSimPlugin.cmake`, after line 112 add:
```cmake
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/MTRSimDriver.cpp
```
and after line 126 add:
```cmake
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/MTRSimDriver.hpp
```

- [ ] **Step 3: Write the failing test**

Create `tests/test_mtrsim_driver.cpp`:

```cpp
#include "LibMTRSim/MTRSimDriver.hpp"

#include <catch2/catch.hpp>

#include <numbers>

TEST_CASE("buildUniformODF produces correct bin centres", "[mtrsim_driver]") {
  const mtrsim::ODFComponent uni = mtrsim::buildUniformODF(72, 36, 72);

  REQUIRE(uni.odfVal.size() == 72 * 36 * 72);
  REQUIRE(uni.phi1Bins.size() == 72 * 36 * 72);

  // Uniform mass: every bin equal, sums to 1.
  REQUIRE(uni.odfVal.sum() == Catch::Approx(1.0));
  REQUIRE(uni.odfVal[0] == Catch::Approx(1.0 / (72.0 * 36.0 * 72.0)));

  // First bin centre: i1=iPHI=i2=0 -> all 0.5 * step.
  REQUIRE(uni.phi1Bins[0] == Catch::Approx(0.5 * 2.0 * std::numbers::pi / 72.0));
  REQUIRE(uni.phiBins[0] == Catch::Approx(0.5 * std::numbers::pi / 36.0));
  REQUIRE(uni.phi2Bins[0] == Catch::Approx(0.5 * 2.0 * std::numbers::pi / 72.0));
}
```

Add `test_mtrsim_driver.cpp` to the test sources in `tests/CMakeLists.txt`
(mirror how `test_odf_sampler.cpp` is listed).

- [ ] **Step 4: Run to verify it fails**

Run: `ctest --test-dir <standalone-build> -R mtrsim_driver --output-on-failure`
Expected: build succeeds, test FAILS only if logic is wrong — if it passes
immediately, good (this helper is straightforward). If the target does not yet
exist, configure first: `cmake --preset <local-standalone-preset>`.

- [ ] **Step 5: Run to verify it passes**

Run: `ctest --test-dir <standalone-build> -R mtrsim_driver --output-on-failure`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/LibMTRSim/MTRSimDriver.hpp src/LibMTRSim/MTRSimDriver.cpp \
        LibMTRSim.cmake MTRSimPlugin.cmake \
        tests/test_mtrsim_driver.cpp tests/CMakeLists.txt
git commit -m "feat(lib): add buildUniformODF grid helper to MTRSimDriver"
```

---

### Task 2: `gridToODFComponent()` — reconstruct an ODF component from grid arrays

**Files:**
- Modify: `src/LibMTRSim/MTRSimDriver.hpp`, `src/LibMTRSim/MTRSimDriver.cpp`
- Modify: `tests/test_mtrsim_driver.cpp`

- [ ] **Step 1: Declare the helper in the header**

Add to `MTRSimDriver.hpp` inside `namespace mtrsim`:

```cpp
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
```

- [ ] **Step 2: Write the failing test**

Add to `tests/test_mtrsim_driver.cpp`:

```cpp
TEST_CASE("gridToODFComponent derives bin centres in radians and normalizes", "[mtrsim_driver]") {
  const int n1 = 72, nPHI = 36, n2 = 72;
  std::vector<double> values(n1 * nPHI * n2, 2.0); // unnormalized constant

  const mtrsim::ODFComponent c =
      mtrsim::gridToODFComponent(values, n1, nPHI, n2, 5.0, 5.0, 5.0);

  REQUIRE(c.odfVal.size() == n1 * nPHI * n2);
  REQUIRE(c.odfVal.sum() == Catch::Approx(1.0)); // normalized
  // 5 deg step -> first bin centre 2.5 deg in radians.
  const double deg2rad = std::numbers::pi / 180.0;
  REQUIRE(c.phi1Bins[0] == Catch::Approx(2.5 * deg2rad));
  REQUIRE(c.phiBins[0] == Catch::Approx(2.5 * deg2rad));
  REQUIRE(c.phi2Bins[0] == Catch::Approx(2.5 * deg2rad));
}
```

- [ ] **Step 3: Run to verify it fails**

Run: `ctest --test-dir <standalone-build> -R mtrsim_driver --output-on-failure`
Expected: FAIL to compile/link ("undefined reference to gridToODFComponent").

- [ ] **Step 4: Implement**

Add to `MTRSimDriver.cpp` (and `#include <numbers>` is already present):

```cpp
ODFComponent gridToODFComponent(const std::vector<double>& values, int n1, int nPHI, int n2, double stepDeg1, double stepDegPHI, double stepDeg2) {
  const int nTotal = n1 * nPHI * n2;
  const double deg2rad = std::numbers::pi / 180.0;
  const double s1 = stepDeg1 * deg2rad;
  const double sP = stepDegPHI * deg2rad;
  const double s2 = stepDeg2 * deg2rad;

  Eigen::VectorXd phi1Bins(nTotal);
  Eigen::VectorXd phiBins(nTotal);
  Eigen::VectorXd phi2Bins(nTotal);
  for (int ix = 0; ix < nTotal; ++ix) {
    const int i1 = ix / (nPHI * n2);
    const int iPHI = (ix % (nPHI * n2)) / n2;
    const int i2 = ix % n2;
    phi1Bins[ix] = (i1 + 0.5) * s1;
    phiBins[ix] = (iPHI + 0.5) * sP;
    phi2Bins[ix] = (i2 + 0.5) * s2;
  }

  ODFComponent c;
  c.odfVal = Eigen::Map<const Eigen::VectorXd>(values.data(), nTotal);
  const double total = c.odfVal.sum();
  if (total > 0.0) {
    c.odfVal /= total;
  }
  c.phi1Bins = std::move(phi1Bins);
  c.phiBins = std::move(phiBins);
  c.phi2Bins = std::move(phi2Bins);
  return c;
}
```

- [ ] **Step 5: Run to verify it passes**

Run: `ctest --test-dir <standalone-build> -R mtrsim_driver --output-on-failure`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/LibMTRSim/MTRSimDriver.hpp src/LibMTRSim/MTRSimDriver.cpp tests/test_mtrsim_driver.cpp
git commit -m "feat(lib): add gridToODFComponent reconstruction helper"
```

---

### Task 3: `remapSimToZYX()` — convert sim voxel order to SIMPLNX order

**Files:**
- Modify: `src/LibMTRSim/MTRSimDriver.hpp`, `src/LibMTRSim/MTRSimDriver.cpp`
- Modify: `tests/test_mtrsim_driver.cpp`

- [ ] **Step 1: Declare the helper**

Add to `MTRSimDriver.hpp`:

```cpp
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
std::vector<T> remapSimToZYX(const std::vector<T>& in, int nx, int ny, int nz) {
  std::vector<T> out(in.size());
  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {
        const std::size_t kSim = static_cast<std::size_t>(iz) * nx * ny + static_cast<std::size_t>(ix) * ny + iy;
        const std::size_t kNx = static_cast<std::size_t>(iz) * ny * nx + static_cast<std::size_t>(iy) * nx + ix;
        out[kNx] = in[kSim];
      }
    }
  }
  return out;
}
```

> Template lives in the header (no .cpp entry needed). Add `#include <cstddef>`
> and `#include <vector>` if not already present.

- [ ] **Step 2: Write the failing test**

Add to `tests/test_mtrsim_driver.cpp`:

```cpp
TEST_CASE("remapSimToZYX moves y-fastest data to x-fastest layout", "[mtrsim_driver]") {
  // 2x3x2 grid (nx=2, ny=3, nz=2). Fill sim-order vector with its own index.
  const int nx = 2, ny = 3, nz = 2;
  std::vector<int> in(nx * ny * nz);
  for (int i = 0; i < nx * ny * nz; ++i) { in[i] = i; }

  const std::vector<int> out = mtrsim::remapSimToZYX(in, nx, ny, nz);

  // Spot-check: SIMPLNX (ix=1, iy=0, iz=0) -> kNx = 1.
  //   source sim index = iz*(nx*ny) + ix*ny + iy = 0 + 1*3 + 0 = 3.
  REQUIRE(out[1] == 3);
  // SIMPLNX (ix=0, iy=1, iz=0) -> kNx = 2; sim = 0 + 0 + 1 = 1.
  REQUIRE(out[2] == 1);
  // Same total size, same multiset of values.
  REQUIRE(out.size() == in.size());
}
```

- [ ] **Step 3: Run to verify it fails, then passes**

Run: `ctest --test-dir <standalone-build> -R mtrsim_driver --output-on-failure`
Expected: FAIL before adding the template (compile error), PASS after.

- [ ] **Step 4: Commit**

```bash
git add src/LibMTRSim/MTRSimDriver.hpp tests/test_mtrsim_driver.cpp
git commit -m "feat(lib): add remapSimToZYX voxel-order helper"
```

---

### Task 4: `simulateMTR()` end-to-end driver (returns SIMPLNX-ordered results)

**Files:**
- Modify: `src/LibMTRSim/MTRSimDriver.hpp`, `src/LibMTRSim/MTRSimDriver.cpp`
- Modify: `src/app/main.cpp`
- Modify: `tests/test_mtrsim_driver.cpp`

- [ ] **Step 1: Declare `MTRSimResult` and `simulateMTR` in the header**

Add to `MTRSimDriver.hpp`:

```cpp
/**
 * @brief Per-voxel simulation output in SIMPLNX z,y,x order.
 */
struct LIBMTRSIM_EXPORT MTRSimResult {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  std::vector<int32_t> mtrIndex;   ///< 1-based component id per voxel, length N
  std::vector<double> phi1;        ///< Euler phi1 [rad] per voxel, length N
  std::vector<double> phi;         ///< Euler PHI  [rad] per voxel, length N
  std::vector<double> phi2;        ///< Euler phi2 [rad] per voxel, length N
};

/**
 * @brief Run the full MTR simulation: PGRF assignment -> per-component ODF
 *        sampling -> per-voxel orientation assignment, returned in SIMPLNX
 *        z,y,x voxel order.
 *
 * @param params           Fully populated SimulationParams (sizes/spacing in a
 *                         single consistent length unit; volumeFractions define
 *                         numComponents; thetaList has >= numComponents-1 rows).
 * @param odfComponents    One ODFComponent per volume-fraction entry, on a
 *                         shared (n1 x nPHI x n2) grid.
 * @param rng              Seeded RNG (mt19937_64).
 * @param n1,nPHI,n2       Bin counts of the ODF grid (for the uniform reference).
 */
LIBMTRSIM_EXPORT MTRSimResult simulateMTR(const SimulationParams& params, const std::vector<ODFComponent>& odfComponents, std::mt19937_64& rng, int n1, int nPHI, int n2);
```

- [ ] **Step 2: Implement `simulateMTR` (logic lifted from `main.cpp:231-319`)**

Add to `MTRSimDriver.cpp` the includes and function:

```cpp
#include "PGRFSimulation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
```

```cpp
MTRSimResult simulateMTR(const SimulationParams& params, const std::vector<ODFComponent>& odfComponents, std::mt19937_64& rng, int n1, int nPHI, int n2) {
  const int nx = static_cast<int>(std::round(params.xLen / params.dx));
  const int ny = static_cast<int>(std::round(params.yLen / params.dy));
  const int nz = std::max(static_cast<int>(std::round(params.zLen / params.dz)), 1);
  const int N = nx * ny * nz;

  if (static_cast<int>(odfComponents.size()) != static_cast<int>(params.volumeFractions.size())) {
    throw std::invalid_argument("simulateMTR: odfComponents count must equal volumeFractions count");
  }

  // 1. PGRF assignment (sim-ordered, 1-based component ids).
  PGRFSimulation pgrf{rng};
  const PGRFResult pgrf_result = pgrf.run(params); // throws on bad dims

  // 2. Sample N orientations per component against the uniform reference.
  const ODFComponent uniformOdf = buildUniformODF(n1, nPHI, n2);
  const int numComponents = static_cast<int>(odfComponents.size());
  std::vector<Eigen::MatrixXd> orientSamples(static_cast<std::size_t>(numComponents));
  ODFSampler sampler{rng};
  for (int j = 0; j < numComponents; ++j) {
    orientSamples[static_cast<std::size_t>(j)] = sampler.sampleN(N, odfComponents[static_cast<std::size_t>(j)], uniformOdf);
  }

  // 3. Assign per-voxel orientation by component (sim order).
  std::vector<int32_t> mtrSim(N);
  std::vector<double> phi1Sim(N), phiSim(N), phi2Sim(N);
  for (int i = 0; i < N; ++i) {
    const int comp = pgrf_result.mtrIndex[i] - 1;
    mtrSim[i] = pgrf_result.mtrIndex[i];
    phi1Sim[i] = orientSamples[static_cast<std::size_t>(comp)](i, 0);
    phiSim[i] = orientSamples[static_cast<std::size_t>(comp)](i, 1);
    phi2Sim[i] = orientSamples[static_cast<std::size_t>(comp)](i, 2);
  }

  // 4. Remap to SIMPLNX z,y,x order.
  MTRSimResult out;
  out.nx = nx; out.ny = ny; out.nz = nz;
  out.mtrIndex = remapSimToZYX(mtrSim, nx, ny, nz);
  out.phi1 = remapSimToZYX(phi1Sim, nx, ny, nz);
  out.phi = remapSimToZYX(phiSim, nx, ny, nz);
  out.phi2 = remapSimToZYX(phi2Sim, nx, ny, nz);
  return out;
}
```

- [ ] **Step 3: Refactor `src/app/main.cpp` to call `simulateMTR`**

Replace the inline blocks `main.cpp:263-319` (PGRF run, ODF load already above,
sampling, per-voxel assign) with a single call. Keep the existing
`loadODFComponents()` and the CSV/PNG output. After loading `odfComponents`,
replace the sampling/assignment with:

```cpp
  mtrsim::MTRSimResult sim = mtrsim::simulateMTR(params, odfComponents, rng, 72, 36, 72);
  // main.cpp writes its CSV/PNG in MATLAB (z,x,y) order historically; for the
  // standalone tool, rebuild phi vectors from sim (now z,y,x) — acceptable, the
  // CSV is diagnostic only. Map sim.phi1/phi/phi2 into the existing Eigen
  // VectorXd phi1Vec/phiVec/phi2Vec used by the writers.
  Eigen::VectorXd phi1Vec = Eigen::Map<Eigen::VectorXd>(sim.phi1.data(), sim.phi1.size());
  Eigen::VectorXd phiVec = Eigen::Map<Eigen::VectorXd>(sim.phi.data(), sim.phi.size());
  Eigen::VectorXd phi2Vec = Eigen::Map<Eigen::VectorXd>(sim.phi2.data(), sim.phi2.size());
```

Remove the now-unused inline `buildUniformODF` from `main.cpp` (it lives in
`MTRSimDriver` now); keep `loadODFComponents`. Update `result.mtrIndex[i]` CSV
references to `sim.mtrIndex[i]`.

> The standalone CSV/PNG are diagnostic; their exact voxel order does not affect
> the filter. The LibMTRSim statistical tests (Step 4) are the source of truth.

- [ ] **Step 4: Write the statistical end-to-end test**

Add to `tests/test_mtrsim_driver.cpp`. Reuse the committed ODF exemplar at
`data/simulation_ODF.h5` via the existing `mtrsim::readODFComponents` (from
`ODFFileIO.hpp`) OR build synthetic components if the test harness lacks a path;
use the project test-data macro already used by `test_odf_sampler.cpp` to locate
`data/`.

```cpp
#include "LibMTRSim/ODFFileIO.hpp"   // mtrsim::readODFComponents (if used)

TEST_CASE("simulateMTR reproduces target volume fractions (statistical)", "[mtrsim_driver][statistical]") {
  mtrsim::SimulationParams params;
  params.xLen = 2.0; params.yLen = 2.0; params.zLen = 0.0;
  params.dx = 0.02; params.dy = 0.02; params.dz = 0.02;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.45, 0.10}, {0.08, 0.37, 0.08}};
  params.seed = 42;

  // Three identical uniform components (sampling correctness is covered by
  // test_odf_sampler.cpp; here we validate assignment fractions + ordering).
  std::vector<mtrsim::ODFComponent> comps = {
      mtrsim::buildUniformODF(72, 36, 72),
      mtrsim::buildUniformODF(72, 36, 72),
      mtrsim::buildUniformODF(72, 36, 72)};

  std::mt19937_64 rng(params.seed);
  const mtrsim::MTRSimResult r = mtrsim::simulateMTR(params, comps, rng, 72, 36, 72);

  const int N = r.nx * r.ny * r.nz;
  REQUIRE(static_cast<int>(r.mtrIndex.size()) == N);

  // MTR ids are exactly the set {1,2,3}.
  for (int v : r.mtrIndex) { REQUIRE(v >= 1); REQUIRE(v <= 3); }

  // Empirical volume fractions within tolerance of targets.
  std::array<int, 3> counts{0, 0, 0};
  for (int v : r.mtrIndex) { counts[v - 1]++; }
  REQUIRE(static_cast<double>(counts[0]) / N == Catch::Approx(0.30).margin(0.05));
  REQUIRE(static_cast<double>(counts[1]) / N == Catch::Approx(0.35).margin(0.05));
  REQUIRE(static_cast<double>(counts[2]) / N == Catch::Approx(0.35).margin(0.05));

  // Euler ranges valid.
  for (double a : r.phi1) { REQUIRE(a >= 0.0); REQUIRE(a <= 2.0 * std::numbers::pi); }
  for (double a : r.phi) { REQUIRE(a >= 0.0); REQUIRE(a <= std::numbers::pi); }
}
```

Add `#include <array>` to the test file.

- [ ] **Step 5: Run to verify it fails, then passes**

Run: `ctest --test-dir <standalone-build> -R mtrsim_driver --output-on-failure`
Expected: FAIL until `simulateMTR` is implemented and `main.cpp` compiles; then
PASS. If volume-fraction margins are too tight on a small grid, increase the
domain (e.g. `xLen=yLen=4.0`) rather than loosening past 0.05.

- [ ] **Step 6: Commit**

```bash
git add src/LibMTRSim/MTRSimDriver.hpp src/LibMTRSim/MTRSimDriver.cpp \
        src/app/main.cpp tests/test_mtrsim_driver.cpp
git commit -m "feat(lib): add simulateMTR end-to-end driver; main.cpp uses it"
```

---

## Phase B — The `MTRSimFilter`

### Task 5: Scaffold the algorithm (`MTRSim`) header + input struct

**Files:**
- Create: `src/MTRSim/Filters/Algorithms/MTRSim.hpp`
- Create: `src/MTRSim/Filters/Algorithms/MTRSim.cpp`

- [ ] **Step 1: Create the algorithm header**

Create `src/MTRSim/Filters/Algorithms/MTRSim.hpp` (mirror
`Algorithms/ReadMTRSimODF.hpp`):

```cpp
#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include <vector>

namespace nx::core
{

struct MTRSIM_EXPORT MTRSimInputValues
{
  DataPath inputOdfGeometryPath;
  std::vector<DataPath> odfComponentPaths;
  std::vector<std::vector<double>> volumeFractions; // 1 row x N cols
  std::vector<std::vector<double>> thetaList;       // M rows x 3 cols
  std::vector<float> physicalSize;                  // [x,y,z] microns
  std::vector<float> physicalSpacing;               // [x,y,z] microns
  uint64 seed;
  bool generatePolarColoring;
  DataPath outputGeometryPath;
  std::string cellAttrMatName;
  std::string mtrIdsArrayName;
  std::string eulersArrayName;
  std::string polarColorsArrayName;
};

/**
 * @class MTRSim
 * @brief Runs the MTR plurigaussian simulation and writes MTR ids, Euler
 * angles, and optional polar-color cell data into the output Image Geometry.
 */
class MTRSIM_EXPORT MTRSim
{
public:
  MTRSim(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, MTRSimInputValues* inputValues);
  ~MTRSim() noexcept;

  MTRSim(const MTRSim&) = delete;
  MTRSim(MTRSim&&) noexcept = delete;
  MTRSim& operator=(const MTRSim&) = delete;
  MTRSim& operator=(MTRSim&&) noexcept = delete;

  Result<> operator()();

private:
  DataStructure& m_DataStructure;
  const MTRSimInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
```

- [ ] **Step 2: Create the algorithm .cpp stub (compiles, returns ok)**

Create `src/MTRSim/Filters/Algorithms/MTRSim.cpp` with the constructor/destructor
and a minimal `operator()` returning `{}` (full body added in Task 7). Mirror
`Algorithms/ReadMTRSimODF.cpp` constructor wiring.

- [ ] **Step 3: Commit**

```bash
git add src/MTRSim/Filters/Algorithms/MTRSim.hpp src/MTRSim/Filters/Algorithms/MTRSim.cpp
git commit -m "feat(filter): scaffold MTRSim algorithm header + stub"
```

---

### Task 6: Scaffold the filter (`MTRSimFilter`) with parameters + preflight

**Files:**
- Create: `src/MTRSim/Filters/MTRSimFilter.hpp`
- Create: `src/MTRSim/Filters/MTRSimFilter.cpp`
- Modify: `MTRSimPlugin.cmake` (`FilterList`, `AlgorithmList`)
- Create/Modify: `test/MTRSimTest.cpp`, `test/CMakeLists.txt`

- [ ] **Step 1: Create the filter header**

Create `src/MTRSim/Filters/MTRSimFilter.hpp` mirroring `ReadMTRSimODFFilter.hpp`,
with these parameter keys and a fresh UUID (generate with `uuidgen`):

```cpp
  static inline constexpr StringLiteral k_InputOdfGeometry_Key = "input_odf_geometry";
  static inline constexpr StringLiteral k_OdfComponentArrays_Key = "odf_component_arrays";
  static inline constexpr StringLiteral k_VolumeFractions_Key = "volume_fractions";
  static inline constexpr StringLiteral k_ThetaList_Key = "theta_list";
  static inline constexpr StringLiteral k_PhysicalSize_Key = "physical_size";
  static inline constexpr StringLiteral k_PhysicalSpacing_Key = "physical_spacing";
  static inline constexpr StringLiteral k_UseSeed_Key = "use_seed";
  static inline constexpr StringLiteral k_SeedValue_Key = "seed_value";
  static inline constexpr StringLiteral k_SeedArrayName_Key = "seed_array_name";
  static inline constexpr StringLiteral k_GeneratePolarColoring_Key = "generate_polar_coloring";
  static inline constexpr StringLiteral k_OutputGeometry_Key = "output_geometry";
  static inline constexpr StringLiteral k_CellAttrMatName_Key = "cell_attribute_matrix_name";
  static inline constexpr StringLiteral k_MtrIdsArrayName_Key = "mtr_ids_array_name";
  static inline constexpr StringLiteral k_EulersArrayName_Key = "eulers_array_name";
  static inline constexpr StringLiteral k_PolarColorsArrayName_Key = "polar_colors_array_name";
```

End the header with:
```cpp
SIMPLNX_DEF_FILTER_TRAITS(nx::core, MTRSimFilter, "<NEW-UUID-FROM-uuidgen>");
```

- [ ] **Step 2: Implement `parameters()`**

In `MTRSimFilter.cpp`, includes (add to the `ReadMTRSimODFFilter.cpp` set):

```cpp
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/MultiArraySelectionParameter.hpp"
#include "simplnx/Parameters/DynamicTableParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Utilities/SIMPLConversion.hpp"
#include "simplnx/Common/Types.hpp"
#include <random>
```

`parameters()` body:

```cpp
Parameters MTRSimFilter::parameters() const
{
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input ODF"});
  params.insert(std::make_unique<GeometrySelectionParameter>(k_InputOdfGeometry_Key, "Input ODF Geometry", "Image Geometry holding the ODF (from the Read/Compute ODF filters).",
                                                             DataPath{}, GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<MultiArraySelectionParameter>(k_OdfComponentArrays_Key, "ODF Component Arrays", "Ordered list of per-component ODF cell arrays. Order maps to Volume Fraction columns.",
                                                              MultiArraySelectionParameter::ValueType{}, MultiArraySelectionParameter::AllowedTypes{IArray::ArrayType::DataArray},
                                                              MultiArraySelectionParameter::AllowedComponentShapes{{1}}, GetAllNumericTypes()));

  params.insertSeparator(Parameters::Separator{"Simulation Parameters"});
  {
    DynamicTableInfo vfInfo;
    vfInfo.setRowsInfo(DynamicTableInfo::StaticVectorInfo(1));
    vfInfo.setColsInfo(DynamicTableInfo::DynamicVectorInfo(1, 3, "Comp {}"));
    params.insert(std::make_unique<DynamicTableParameter>(k_VolumeFractions_Key, "Volume Fraction", "One value per ODF component; must match the component count and sum to 1.0.", vfInfo));
  }
  {
    DynamicTableInfo thetaInfo;
    thetaInfo.setRowsInfo(DynamicTableInfo::DynamicVectorInfo(1, 2, "Gaussian {}"));
    thetaInfo.setColsInfo(DynamicTableInfo::StaticVectorInfo({"theta_x", "theta_y", "theta_z"}));
    params.insert(std::make_unique<DynamicTableParameter>(k_ThetaList_Key, "Theta List", "Correlation lengths [theta_x, theta_y, theta_z] per latent Gaussian. Needs >= (components - 1) rows. Same length unit as Physical Size/Spacing.", thetaInfo));
  }
  params.insert(std::make_unique<VectorFloat32Parameter>(k_PhysicalSize_Key, "Physical Size (microns)", "Domain extent X,Y,Z.", std::vector<float32>{38.1f, 12.7f, 0.0f}, std::vector<std::string>{"X", "Y", "Z"}));
  params.insert(std::make_unique<VectorFloat32Parameter>(k_PhysicalSpacing_Key, "Physical Spacing (microns)", "Voxel spacing X,Y,Z.", std::vector<float32>{0.02f, 0.02f, 0.02f}, std::vector<std::string>{"X", "Y", "Z"}));

  params.insertSeparator(Parameters::Separator{"Random Number Seed Parameters"});
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_UseSeed_Key, "Use Seed for Random Generation", "When true the user can supply a fixed seed.", false));
  params.insert(std::make_unique<NumberParameter<uint64>>(k_SeedValue_Key, "Seed Value", "The seed fed into the random generator.", std::mt19937::default_seed));
  params.insert(std::make_unique<DataObjectNameParameter>(k_SeedArrayName_Key, "Stored Seed Value Array Name", "Top-level array recording the seed used.", "MTRSim SeedValue"));

  params.insertSeparator(Parameters::Separator{"Outputs"});
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_GeneratePolarColoring_Key, "Generate Polar Coloring", "Create a 3-component UInt8 RGB array using the MATLAB polar color mapping.", false));
  params.insert(std::make_unique<DataGroupCreationParameter>(k_OutputGeometry_Key, "Output Image Geometry", "Path of the new microstructure Image Geometry.", DataPath({"MTR Microstructure"})));
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellAttrMatName_Key, "Cell Attribute Matrix Name", "Name of the created cell AttributeMatrix.", "Cell Data"));
  params.insert(std::make_unique<DataObjectNameParameter>(k_MtrIdsArrayName_Key, "MTR Ids Array Name", "Int32 per-voxel MTR component id (1-based).", "MTRIds"));
  params.insert(std::make_unique<DataObjectNameParameter>(k_EulersArrayName_Key, "Euler Angles Array Name", "Float32 3-component Bunge Euler angles [rad].", "Eulers"));
  params.insert(std::make_unique<DataObjectNameParameter>(k_PolarColorsArrayName_Key, "Polar Colors Array Name", "UInt8 3-component RGB polar coloring.", "Polar Colors"));

  params.linkParameters(k_UseSeed_Key, k_SeedValue_Key, true);
  params.linkParameters(k_GeneratePolarColoring_Key, k_PolarColorsArrayName_Key, true);

  return params;
}
```

> Add `#include "simplnx/Parameters/util/DynamicTableInfo.hpp"`. Verify
> `GetAllNumericTypes()`/`VectorFloat32Parameter` include paths against
> `ReadMTRSimODFFilter.cpp` siblings; adjust to the names the local simplnx
> exposes (`GetAllNumericTypes` is in `MultiArraySelectionParameter.hpp` usage
> in other filters — grep `simplnx` if the symbol differs).

- [ ] **Step 3: Implement `preflightImpl` (validation + create geometry/arrays)**

```cpp
IFilter::PreflightResult MTRSimFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                     const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  auto pOdfArrays = filterArgs.value<MultiArraySelectionParameter::ValueType>(k_OdfComponentArrays_Key);
  auto pVolumeFractions = filterArgs.value<DynamicTableParameter::ValueType>(k_VolumeFractions_Key);
  auto pThetaList = filterArgs.value<DynamicTableParameter::ValueType>(k_ThetaList_Key);
  auto pSize = filterArgs.value<std::vector<float32>>(k_PhysicalSize_Key);
  auto pSpacing = filterArgs.value<std::vector<float32>>(k_PhysicalSpacing_Key);
  auto pUseSeed = filterArgs.value<bool>(k_UseSeed_Key);
  auto pGenPolar = filterArgs.value<bool>(k_GeneratePolarColoring_Key);
  auto pOutGeomPath = filterArgs.value<DataPath>(k_OutputGeometry_Key);
  auto pCellAttrMatName = filterArgs.value<std::string>(k_CellAttrMatName_Key);
  auto pMtrIdsName = filterArgs.value<std::string>(k_MtrIdsArrayName_Key);
  auto pEulersName = filterArgs.value<std::string>(k_EulersArrayName_Key);
  auto pPolarName = filterArgs.value<std::string>(k_PolarColorsArrayName_Key);
  auto pSeedArrayName = filterArgs.value<std::string>(k_SeedArrayName_Key);

  nx::core::Result<OutputActions> resultOutputActions;
  std::vector<PreflightValue> preflightUpdatedValues;

  const usize numComponents = pOdfArrays.size();
  if (numComponents < 2)
  {
    return {MakeErrorResult<OutputActions>(-13001, "MTRSim requires at least 2 ODF component arrays.")};
  }

  // Volume fractions: exactly 1 row, numComponents columns, sum ~ 1.0.
  if (pVolumeFractions.size() != 1 || pVolumeFractions[0].size() != numComponents)
  {
    return {MakeErrorResult<OutputActions>(-13002, fmt::format("Volume Fraction must be 1 row x {} columns (one per ODF component).", numComponents))};
  }
  double vfSum = 0.0;
  for (double v : pVolumeFractions[0]) { vfSum += v; }
  if (std::abs(vfSum - 1.0) > 1.0e-3)
  {
    return {MakeErrorResult<OutputActions>(-13003, fmt::format("Volume Fraction values must sum to 1.0 (got {:.4f}).", vfSum))};
  }

  // Theta list: >= numComponents - 1 rows, 3 columns each.
  if (pThetaList.size() < numComponents - 1)
  {
    return {MakeErrorResult<OutputActions>(-13004, fmt::format("Theta List needs at least {} rows (components - 1).", numComponents - 1))};
  }
  for (const auto& row : pThetaList)
  {
    if (row.size() != 3)
    {
      return {MakeErrorResult<OutputActions>(-13005, "Each Theta List row must have exactly 3 columns.")};
    }
  }

  // Grid dims from size/spacing.
  const auto dim = [](float len, float sp) { return static_cast<usize>(std::max(std::lround(len / sp), 1L)); };
  const usize nx = dim(pSize[0], pSpacing[0]);
  const usize ny = dim(pSize[1], pSpacing[1]);
  const usize nz = (pSize[2] <= 0.0f) ? 1 : dim(pSize[2], pSpacing[2]);

  const std::vector<usize> imageGeomDimsXYZ = {nx, ny, nz};
  const std::vector<float32> origin = {0.0f, 0.0f, 0.0f};
  const std::vector<float32> spacingXYZ = {pSpacing[0], pSpacing[1], pSpacing[2]};
  const std::vector<usize> tupleShapeZYX = {nz, ny, nx};

  resultOutputActions.value().appendAction(std::make_unique<CreateImageGeometryAction>(pOutGeomPath, imageGeomDimsXYZ, origin, spacingXYZ, pCellAttrMatName));

  const DataPath cellAttrMatPath = pOutGeomPath.createChildPath(pCellAttrMatName);
  resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::int32, tupleShapeZYX, std::vector<usize>{1}, cellAttrMatPath.createChildPath(pMtrIdsName)));
  resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::float32, tupleShapeZYX, std::vector<usize>{3}, cellAttrMatPath.createChildPath(pEulersName)));
  if (pGenPolar)
  {
    resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::uint8, tupleShapeZYX, std::vector<usize>{3}, cellAttrMatPath.createChildPath(pPolarName)));
  }

  // Seed array (top-level UInt64 scalar).
  resultOutputActions.value().appendAction(std::make_unique<CreateArrayAction>(DataType::uint64, std::vector<usize>{1}, std::vector<usize>{1}, DataPath({pSeedArrayName})));

  preflightUpdatedValues.push_back({"Output Grid (X, Y, Z)", fmt::format("{} x {} x {}", nx, ny, nz)});
  preflightUpdatedValues.push_back({"Number of ODF Components", std::to_string(numComponents)});

  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}
```

> Includes needed: `CreateImageGeometryAction.hpp`, `CreateArrayAction.hpp`,
> `<cmath>`, `fmt/format.h` (already in sibling). The
> `name()/className()/uuid()/humanName()/defaultTags()/parametersVersion()/clone()/FromSIMPLJson()`
> methods mirror `ReadMTRSimODFFilter.cpp` exactly. `humanName()` returns
> `"Generate Synthetic Microtexture"`; `defaultTags()` returns
> `{className(), "MTRSim", "Synthetic", "Microtexture", "Generate"}`.

- [ ] **Step 4: Implement `executeImpl` (delegate to algorithm)**

```cpp
Result<> MTRSimFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                   const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  MTRSimInputValues inputValues;
  inputValues.inputOdfGeometryPath = filterArgs.value<DataPath>(k_InputOdfGeometry_Key);
  inputValues.odfComponentPaths = filterArgs.value<MultiArraySelectionParameter::ValueType>(k_OdfComponentArrays_Key);
  inputValues.volumeFractions = filterArgs.value<DynamicTableParameter::ValueType>(k_VolumeFractions_Key);
  inputValues.thetaList = filterArgs.value<DynamicTableParameter::ValueType>(k_ThetaList_Key);
  inputValues.physicalSize = filterArgs.value<std::vector<float32>>(k_PhysicalSize_Key);
  inputValues.physicalSpacing = filterArgs.value<std::vector<float32>>(k_PhysicalSpacing_Key);
  inputValues.generatePolarColoring = filterArgs.value<bool>(k_GeneratePolarColoring_Key);
  inputValues.outputGeometryPath = filterArgs.value<DataPath>(k_OutputGeometry_Key);
  inputValues.cellAttrMatName = filterArgs.value<std::string>(k_CellAttrMatName_Key);
  inputValues.mtrIdsArrayName = filterArgs.value<std::string>(k_MtrIdsArrayName_Key);
  inputValues.eulersArrayName = filterArgs.value<std::string>(k_EulersArrayName_Key);
  inputValues.polarColorsArrayName = filterArgs.value<std::string>(k_PolarColorsArrayName_Key);

  // Standard simplnx seed handling.
  uint64 seed = filterArgs.value<uint64>(k_SeedValue_Key);
  if (!filterArgs.value<bool>(k_UseSeed_Key))
  {
    seed = static_cast<uint64>(std::chrono::steady_clock::now().time_since_epoch().count());
  }
  dataStructure.getDataRefAs<UInt64Array>(DataPath({filterArgs.value<std::string>(k_SeedArrayName_Key)}))[0] = seed;
  inputValues.seed = seed;

  return MTRSim(dataStructure, messageHandler, shouldCancel, &inputValues)();
}
```

> Add `#include <chrono>` and `#include "simplnx/DataStructure/DataArray.hpp"`.

- [ ] **Step 5: Register filter + algorithm in `MTRSimPlugin.cmake`**

In `FilterList` add `MTRSimFilter`; in `AlgorithmList` add `MTRSim`.

- [ ] **Step 6: Add the test file to `test/CMakeLists.txt` and write preflight error tests**

Add `MTRSimTest.cpp` to `${PLUGIN_NAME}UnitTest_SRCS`. Create
`test/MTRSimTest.cpp` (mirror `ReadMTRSimODFTest.cpp` includes/structure):

```cpp
#include "MTRSim/Filters/MTRSimFilter.hpp"

#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Parameters/DynamicTableParameter.hpp"
#include "simplnx/Parameters/MultiArraySelectionParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"

#include <catch2/catch.hpp>

using namespace nx::core;

namespace
{
// Build a DataStructure with an ODF ImageGeom (72x36x72) holding `n` Float64
// single-component cell arrays named component_0.. and return the geometry +
// component paths. Helper used by every test below.
DataStructure makeOdfDataStructure(usize numComponents, std::vector<DataPath>& outCompPaths)
{
  DataStructure ds;
  auto* ig = ImageGeom::Create(ds, "ODF");
  ig->setDimensions({72, 36, 72});
  ig->setSpacing({5.0f, 5.0f, 5.0f});
  ig->setOrigin({0.0f, 0.0f, 0.0f});
  auto* am = AttributeMatrix::Create(ds, "Cell Data", {72 * 36 * 72}, ig->getId());
  for (usize c = 0; c < numComponents; ++c)
  {
    const std::string name = fmt::format("component_{}", c);
    auto* arr = UInt64Array::CreateWithStore<Float64DataStore>(ds, name, {72 * 36 * 72}, {1}, am->getId());
    (void)arr;
    outCompPaths.push_back(DataPath({"ODF", "Cell Data", name}));
  }
  return ds;
}
} // namespace

TEST_CASE("MTRSimFilter: preflight rejects mismatched volume fraction count", "[MTRSim]")
{
  std::vector<DataPath> compPaths;
  DataStructure ds = makeOdfDataStructure(3, compPaths);

  MTRSimFilter filter;
  Arguments args;
  args.insertOrAssign(MTRSimFilter::k_InputOdfGeometry_Key, DataPath({"ODF"}));
  args.insertOrAssign(MTRSimFilter::k_OdfComponentArrays_Key, compPaths);
  args.insertOrAssign(MTRSimFilter::k_VolumeFractions_Key, DynamicTableParameter::ValueType{{0.5, 0.5}}); // only 2, need 3
  args.insertOrAssign(MTRSimFilter::k_ThetaList_Key, DynamicTableParameter::ValueType{{0.1, 0.45, 0.1}, {0.08, 0.37, 0.08}});
  args.insertOrAssign(MTRSimFilter::k_PhysicalSize_Key, std::vector<float32>{2.0f, 2.0f, 0.0f});
  args.insertOrAssign(MTRSimFilter::k_PhysicalSpacing_Key, std::vector<float32>{0.02f, 0.02f, 0.02f});
  args.insertOrAssign(MTRSimFilter::k_UseSeed_Key, true);
  args.insertOrAssign(MTRSimFilter::k_SeedValue_Key, static_cast<uint64>(42));
  args.insertOrAssign(MTRSimFilter::k_SeedArrayName_Key, std::string("MTRSim SeedValue"));
  args.insertOrAssign(MTRSimFilter::k_GeneratePolarColoring_Key, false);
  args.insertOrAssign(MTRSimFilter::k_OutputGeometry_Key, DataPath({"MTR Microstructure"}));
  args.insertOrAssign(MTRSimFilter::k_CellAttrMatName_Key, std::string("Cell Data"));
  args.insertOrAssign(MTRSimFilter::k_MtrIdsArrayName_Key, std::string("MTRIds"));
  args.insertOrAssign(MTRSimFilter::k_EulersArrayName_Key, std::string("Eulers"));
  args.insertOrAssign(MTRSimFilter::k_PolarColorsArrayName_Key, std::string("Polar Colors"));

  auto result = filter.preflight(ds, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(result.outputActions);
}
```

> Fix the obvious paste error in the helper: components must be created as
> `Float64Array::CreateWithStore<Float64DataStore>` (the `UInt64Array` token is a
> typo — use `Float64Array`). Reuse this `args`-builder as a lambda in later
> tests to stay DRY.

- [ ] **Step 7: Build and run the preflight test**

Run: `cmake --build /Users/mjackson/Workspace7/simplnx/build --config Release`
then `ctest --test-dir /Users/mjackson/Workspace7/simplnx/build -R MTRSim --output-on-failure`
Expected: the mismatch test PASSES (preflight returns invalid). Also add and run
analogous tests for: theta rows `< numComponents-1` (error -13004), VF sum ≠ 1
(error -13003). Each follows the same builder with one field changed.

- [ ] **Step 8: Commit**

```bash
git add src/MTRSim/Filters/MTRSimFilter.hpp src/MTRSim/Filters/MTRSimFilter.cpp \
        MTRSimPlugin.cmake test/MTRSimTest.cpp test/CMakeLists.txt
git commit -m "feat(filter): add MTRSimFilter params + preflight validation + error tests"
```

---

### Task 7: Algorithm `operator()` — read ODF, run sim, write arrays

**Files:**
- Modify: `src/MTRSim/Filters/Algorithms/MTRSim.cpp`
- Modify: `test/MTRSimTest.cpp`

- [ ] **Step 1: Implement `MTRSim::operator()`**

```cpp
#include "MTRSim.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"

#include "LibMTRSim/MTRSimDriver.hpp"

#include <fmt/format.h>

#include <cmath>
#include <exception>

using namespace nx::core;

MTRSim::MTRSim(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, MTRSimInputValues* inputValues)
: m_DataStructure(dataStructure), m_InputValues(inputValues), m_ShouldCancel(shouldCancel), m_MessageHandler(mesgHandler)
{
}

MTRSim::~MTRSim() noexcept = default;

Result<> MTRSim::operator()()
{
  // 1. Pull ODF bin geometry (degrees) from the input ODF ImageGeom.
  const auto& odfGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->inputOdfGeometryPath);
  const SizeVec3 odfDims = odfGeom.getDimensions();    // X=phi2, Y=PHI, Z=phi1
  const FloatVec3 odfSpacing = odfGeom.getSpacing();   // degrees
  const int n2 = static_cast<int>(odfDims[0]);
  const int nPHI = static_cast<int>(odfDims[1]);
  const int n1 = static_cast<int>(odfDims[2]);

  // 2. Reconstruct ODFComponents from the selected cell arrays.
  std::vector<mtrsim::ODFComponent> components;
  components.reserve(m_InputValues->odfComponentPaths.size());
  for (const auto& path : m_InputValues->odfComponentPaths)
  {
    const auto& arr = m_DataStructure.getDataRefAs<Float64Array>(path);
    const auto& store = arr.getDataStoreRef();
    std::vector<double> values(store.begin(), store.end());
    components.push_back(mtrsim::gridToODFComponent(values, n1, nPHI, n2,
        static_cast<double>(odfSpacing[2]), static_cast<double>(odfSpacing[1]), static_cast<double>(odfSpacing[0])));
  }

  // 3. Build SimulationParams.
  mtrsim::SimulationParams params;
  params.xLen = m_InputValues->physicalSize[0];
  params.yLen = m_InputValues->physicalSize[1];
  params.zLen = m_InputValues->physicalSize[2];
  params.dx = m_InputValues->physicalSpacing[0];
  params.dy = m_InputValues->physicalSpacing[1];
  params.dz = m_InputValues->physicalSpacing[2];
  params.volumeFractions = m_InputValues->volumeFractions[0]; // 1 row
  params.thetaList = m_InputValues->thetaList;
  params.seed = m_InputValues->seed;

  // 4. Run the simulation (returns SIMPLNX z,y,x order).
  std::mt19937_64 rng(m_InputValues->seed);
  mtrsim::MTRSimResult sim;
  try
  {
    sim = mtrsim::simulateMTR(params, components, rng, n1, nPHI, n2);
  } catch (const std::exception& e)
  {
    return MakeErrorResult(-13050, fmt::format("MTR simulation failed: {}", e.what()));
  }

  if (m_ShouldCancel) { return {}; }

  // 5. Write MTR ids + Euler arrays.
  const DataPath cellAm = m_InputValues->outputGeometryPath.createChildPath(m_InputValues->cellAttrMatName);
  auto& mtrIds = m_DataStructure.getDataRefAs<Int32Array>(cellAm.createChildPath(m_InputValues->mtrIdsArrayName));
  auto& eulers = m_DataStructure.getDataRefAs<Float32Array>(cellAm.createChildPath(m_InputValues->eulersArrayName));
  auto& mtrStore = mtrIds.getDataStoreRef();
  auto& eulerStore = eulers.getDataStoreRef();

  const std::size_t N = sim.mtrIndex.size();
  for (std::size_t i = 0; i < N; ++i)
  {
    mtrStore[i] = sim.mtrIndex[i];
    eulerStore[i * 3 + 0] = static_cast<float>(sim.phi1[i]);
    eulerStore[i * 3 + 1] = static_cast<float>(sim.phi[i]);
    eulerStore[i * 3 + 2] = static_cast<float>(sim.phi2[i]);
  }

  // 6. Optional polar coloring (Task 8).
  if (m_InputValues->generatePolarColoring)
  {
    return applyPolarColoring(sim, cellAm); // declared/added in Task 8
  }

  return {};
}
```

> If `applyPolarColoring` is not yet added, temporarily inline a `return {};` for
> the polar branch and replace in Task 8. Confirm DataArray alias names
> (`Int32Array`, `Float32Array`, `Float64Array`) against `simplnx`
> `DataArray.hpp`.

- [ ] **Step 2: Write the end-to-end filter test (statistical + array contract)**

Add to `test/MTRSimTest.cpp`. Reuse the `args`-builder from Task 6 with a valid
3-component VF and a fixed seed; after `execute`, assert array existence, types,
component counts, that ids ∈ {1,2,3}, and empirical volume fractions within
0.06 of targets.

```cpp
TEST_CASE("MTRSimFilter: execute produces valid MTR ids + Eulers", "[MTRSim]")
{
  std::vector<DataPath> compPaths;
  DataStructure ds = makeOdfDataStructure(3, compPaths);
  // Fill each ODF component uniformly so sampling is well-defined.
  for (const auto& p : compPaths)
  {
    auto& a = ds.getDataRefAs<Float64Array>(p);
    a.fill(1.0);
  }

  MTRSimFilter filter;
  Arguments args; /* ... same builder, VF = {{0.30,0.35,0.35}}, UseSeed=true, Seed=42, Size={4,4,0} ... */

  auto pre = filter.preflight(ds, args);
  SIMPLNX_RESULT_REQUIRE_VALID(pre.outputActions);
  auto exec = filter.execute(ds, args);
  SIMPLNX_RESULT_REQUIRE_VALID(exec.result);

  const DataPath cellAm({"MTR Microstructure", "Cell Data"});
  auto& ids = ds.getDataRefAs<Int32Array>(cellAm.createChildPath("MTRIds"));
  auto& eul = ds.getDataRefAs<Float32Array>(cellAm.createChildPath("Eulers"));
  REQUIRE(eul.getNumberOfComponents() == 3);

  std::array<int, 3> counts{0, 0, 0};
  for (usize i = 0; i < ids.getNumberOfTuples(); ++i)
  {
    const int v = ids[i];
    REQUIRE(v >= 1); REQUIRE(v <= 3);
    counts[v - 1]++;
  }
  const double n = static_cast<double>(ids.getNumberOfTuples());
  REQUIRE(counts[0] / n == Catch::Approx(0.30).margin(0.06));
  REQUIRE(counts[1] / n == Catch::Approx(0.35).margin(0.06));

  // Seed array recorded.
  auto& seedArr = ds.getDataRefAs<UInt64Array>(DataPath({"MTRSim SeedValue"}));
  REQUIRE(seedArr[0] == 42);
}
```

- [ ] **Step 3: Build, run, verify**

Run: `cmake --build /Users/mjackson/Workspace7/simplnx/build --config Release`
then `ctest --test-dir /Users/mjackson/Workspace7/simplnx/build -R MTRSim --output-on-failure`
Expected: PASS. If volume-fraction margins fail on a small grid, raise `Size` to
`{6,6,0}`.

- [ ] **Step 4: Commit**

```bash
git add src/MTRSim/Filters/Algorithms/MTRSim.cpp test/MTRSimTest.cpp
git commit -m "feat(filter): implement MTRSim algorithm execute + statistical test"
```

---

### Task 8: Optional polar coloring output

**Files:**
- Modify: `src/MTRSim/Filters/Algorithms/MTRSim.hpp`, `src/MTRSim/Filters/Algorithms/MTRSim.cpp`
- Modify: `test/MTRSimTest.cpp`

- [ ] **Step 1: Declare the helper**

Add to `MTRSim.hpp` private section:

```cpp
  Result<> applyPolarColoring(const struct mtrsim::MTRSimResult& sim, const DataPath& cellAttrMatPath);
```
(Forward declaration alternative: include `LibMTRSim/MTRSimDriver.hpp` in the
header and drop the `struct` keyword. Match whichever compiles cleanly.)

- [ ] **Step 2: Implement using `IPFMapper` (MATLAB/polar scheme)**

Add to `MTRSim.cpp`:

```cpp
#include "LibMTRSim/IPFMapper.hpp"

Result<> MTRSim::applyPolarColoring(const mtrsim::MTRSimResult& sim, const DataPath& cellAttrMatPath)
{
  const std::size_t N = sim.phi1.size();
  Eigen::VectorXd phi1 = Eigen::Map<const Eigen::VectorXd>(sim.phi1.data(), N);
  Eigen::VectorXd phi = Eigen::Map<const Eigen::VectorXd>(sim.phi.data(), N);
  Eigen::VectorXd phi2 = Eigen::Map<const Eigen::VectorXd>(sim.phi2.data(), N);

  mtrsim::IPFMapper mapper{mtrsim::CrystalSystem::HCP};
  const std::vector<mtrsim::RGBColor> colors = mapper.eulerToColors(phi1, phi, phi2, {0.0, 0.0, 1.0}, mtrsim::IPFColorScheme::MatLab);

  auto& rgb = m_DataStructure.getDataRefAs<UInt8Array>(cellAttrMatPath.createChildPath(m_InputValues->polarColorsArrayName));
  auto& store = rgb.getDataStoreRef();
  for (std::size_t i = 0; i < N; ++i)
  {
    store[i * 3 + 0] = colors[i].r;
    store[i * 3 + 1] = colors[i].g;
    store[i * 3 + 2] = colors[i].b;
  }
  return {};
}
```

> `eulerToColors` returns values already in SIMPLNX voxel order because `sim`
> is already remapped. Add `#include <Eigen/Dense>` to the .cpp.

- [ ] **Step 3: Write the polar-coloring test**

Add to `test/MTRSimTest.cpp`: a copy of the execute test with
`k_GeneratePolarColoring_Key = true`, asserting the `Polar Colors` array exists,
is `UInt8`, 3 components, same tuple count as `MTRIds`, and is not all-zero:

```cpp
TEST_CASE("MTRSimFilter: polar coloring produces a populated RGB array", "[MTRSim]")
{
  // ... build ds + args as in the execute test, but:
  args.insertOrAssign(MTRSimFilter::k_GeneratePolarColoring_Key, true);
  // ... preflight VALID, execute VALID ...
  const DataPath cellAm({"MTR Microstructure", "Cell Data"});
  auto& rgb = ds.getDataRefAs<UInt8Array>(cellAm.createChildPath("Polar Colors"));
  REQUIRE(rgb.getNumberOfComponents() == 3);
  REQUIRE(rgb.getNumberOfTuples() > 0);
  uint64 sum = 0;
  for (usize i = 0; i < rgb.getSize(); ++i) { sum += rgb[i]; }
  REQUIRE(sum > 0);
}
```

Also assert that when the bool is **false**, the `Polar Colors` array does NOT
exist:
```cpp
  REQUIRE(ds.getDataAs<UInt8Array>(DataPath({"MTR Microstructure", "Cell Data", "Polar Colors"})) == nullptr);
```

- [ ] **Step 4: Build, run, verify, commit**

Run the build + `ctest -R MTRSim`. Expected PASS.
```bash
git add src/MTRSim/Filters/Algorithms/MTRSim.hpp src/MTRSim/Filters/Algorithms/MTRSim.cpp test/MTRSimTest.cpp
git commit -m "feat(filter): add optional polar coloring output + tests"
```

---

## Phase C — Integration, CI, docs

### Task 9: Full dual-build + run all tests

- [ ] **Step 1: Run clang-format on all new files**

Run the repo's format script (mirror `git log` commit "Clang Format"; the CI
`format_pr.yml` enforces it). Example:
`find src/LibMTRSim src/MTRSim tests test -name 'MTRSim*' -o -name 'test_mtrsim*' | xargs clang-format -i`
(adjust glob). Commit any reformatting.

- [ ] **Step 2: Build + test the standalone library**

Run: `cmake --build <standalone-build>` then
`ctest --test-dir <standalone-build> --output-on-failure`
Expected: all `mtrsim_driver` + existing lib tests PASS.

- [ ] **Step 3: Build + test the plugin**

Run: `cmake --build /Users/mjackson/Workspace7/simplnx/build --config Release`
then `ctest --test-dir /Users/mjackson/Workspace7/simplnx/build -R MTRSim --output-on-failure`
Expected: all filter tests PASS.

- [ ] **Step 4: Commit any fixes**

```bash
git add -A && git commit -m "chore: clang-format + dual-build fixes for MTRSim filter"
```

---

### Task 10: Filter documentation stub

**Files:**
- Create: `docs/MTRSim/MTRSimFilter.md`

- [ ] **Step 1: Write the doc** using the `bluequartz-skills:filter-documentation`
  template: summary, the parameter table (ODF geometry, component arrays, volume
  fraction, theta list, size/spacing, seed group, polar coloring, outputs),
  outputs description (MTRIds 1-based, Eulers radians, Polar Colors), the
  units-consistency note, and an example pipeline reference. (Full polish is
  Milestone AL.)

- [ ] **Step 2: Commit**

```bash
git add docs/MTRSim/MTRSimFilter.md
git commit -m "docs: add MTRSimFilter documentation stub"
```

---

### Task 11: Validate CI

- [ ] **Step 1: Push the branch and open a PR to `develop`**

```bash
git push -u origin topic/create_mtr_sim_filter
gh pr create --base develop --title "Milestone AK: MTRSimFilter (Generate Synthetic Microtexture)" --body "Implements the MTRSim simulation filter, deterministic + statistical tests, and validates CI. See docs/superpowers/plans/2026-06-01-milestone-ak.md."
```

- [ ] **Step 2: Watch the workflows**

Run: `gh pr checks --watch`
Expected: `linux`, `macos`, `windows`, `format_pr` all green. If `format_pr`
fails, run clang-format and push. If a build fails on a missing vcpkg dependency
for the embedded LibMTRSim (Eigen/spdlog/CLI11/nlohmann-json/hdf5/stb), confirm
those appear in `vcpkg.json` and are visible to the simplnx build; add any
missing entry and push.

- [ ] **Step 3: Confirm the acceptance criteria**

Verify the CI run shows the `MTRSim` ctest cases executed and passed on all
platforms (the "all new code paths exercised by passing tests in CI"
requirement). Capture a screenshot/log link for the report.

---

## Self-Review Notes (for the executor)

- **Spec coverage:** filter (Task 5-8), component-exact tests (Task 1-3, 6),
  pipeline-statistical test (Task 4, 7), preflight/error tests (Task 6), CI
  (Task 11), docs stub (Task 10), `buildUniformODF` exposure (Task 1), ODF→
  component helper (Task 2), voxel remap (Task 3), units note (Task 5 param help
  + Task 10 doc). All spec sections map to a task.
- **Known things to verify against the live simplnx headers** (APIs evolve):
  `GetAllNumericTypes`, `VectorFloat32Parameter`, `DynamicTableInfo::DynamicVectorInfo`
  signature, `AttributeMatrix::Create`, `Float64Array::CreateWithStore`,
  `SIMPLNX_RESULT_REQUIRE_VALID/INVALID` macro names. Grep the local checkout at
  `/Users/mjackson/Workspace7/simplnx` and match a sibling filter/test exactly if
  any symbol differs.
- **Fix the deliberate typo** flagged in Task 6 Step 6 (`UInt64Array` →
  `Float64Array` in the ODF test helper).

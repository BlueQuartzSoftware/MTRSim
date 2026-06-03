# Milestone AK — Integration of MTR Representation Codes into DREAM3D-NX Filters

**Date:** 2026-06-01
**Status:** Design approved — ready for implementation plan
**Task reference:** `.claude/mtr_sbir_tasks.md` → Milestone AK (3.2.3 Task 3)

---

## 1. Purpose & Deliverables

Milestone AK exposes the MTR microtexture-representation building capability —
already ported to C++ in `LibMTRSim` (Milestone AH) and wrapped with ODF
import/export/compute filters (Milestone AJ) — to end users through a single
DREAM3D-NX pipeline filter, backed by automated tests and continuous
integration.

Three things ship together:

- **A. The filter** — one new DREAM3D-NX filter that runs the full MTR
  simulation and writes the results into the DataStructure.
- **B. Tests** — unit tests that exercise the *filter's* value-added code paths
  (the library is already tested in `tests/`), including a statistical
  end-to-end test that reuses the LibMTRSim exemplar/seed, plus preflight/error
  tests.
- **C. CI** — GitHub Actions workflows (linux, macOS, Windows, clang-format)
  that build the full `simplnx` + `MTRSim` plugin and run `ctest` automatically
  on push/PR. **Already added** under `.github/`; this milestone validates it
  runs green.

A narrative **report** (per the project-wide convention) is written after the
work is implemented and verified.

### Acceptance criteria mapping

> "All new code paths shall be exercised by passing tests in the CI
> environment, and the framework shall be configured to run automatically upon
> code changes within the private repository."

- New code paths = the `MTRSim` filter + algorithm → covered by **B**.
- "Run automatically upon code changes" → covered by **C** (triggers on
  push/PR to `develop`/`main`).

---

## 2. The Filter

| Field | Value |
|---|---|
| Filter class | `MTRSimFilter` |
| Algorithm class | `MTRSim` (in `Filters/Algorithms/`) |
| Human name | **Generate Synthetic Microtexture** |
| Location | `src/MTRSim/Filters/` |

### 2.1 Parameters (inputs)

| Parameter | Type | Notes / preflight rules |
|---|---|---|
| Input ODF Geometry | `GeometrySelectionParameter` (ImageGeom) | The ODF Image Geometry built by the AJ Read/Compute ODF filters. Bin sizes are read from the geometry's spacing (degrees). |
| ODF Component Arrays | `MultiArraySelectionParameter` (Float64) | Explicit, **ordered** list of the per-component cell arrays. Count defines `numComponents`. Order is significant — index `j` maps to `volumeFractions[j]`. |
| Volume Fraction | `DynamicTableParameter` — 1 fixed row × N cols | Preflight: column count == `numComponents`; values sum ≈ 1.0 (tolerance). Converted to `float`/`double`. |
| Theta List | `DynamicTableParameter` — 3 fixed cols × M rows | Preflight: `M >= numComponents - 1`. Columns are `[theta_x, theta_y, theta_z]` correlation lengths. |
| Physical Size | `VectorFloat32Parameter` (FloatVec3), µm | Domain extent. |
| Physical Spacing | `VectorFloat32Parameter` (FloatVec3), µm | Voxel spacing. |
| Use Seed for Random Generation | `BoolParameter` (default `false`), linkable | Standard simplnx seed pattern (see below). |
| Seed Value | `NumberParameter<uint64>` (default `std::mt19937::default_seed`) | Linked to "Use Seed"; enabled when it is on. |
| Stored Seed Value Array Name | `DataObjectNameParameter` (default `"MTRSim SeedValue"`) | Top-level UInt64 array that records the seed actually used. |
| Generate Polar Coloring | `BoolParameter` (default `false`) | Gates creation of the RGB output array. |
| *nuggetVariance* | — | **Not exposed.** Unused by the simulation. |

> **Random seed pattern.** Follow the established simplnx convention (e.g.
> `MergeTwinsFilter`): a linkable `BoolParameter` "Use Seed for Random
> Generation" gates a `NumberParameter<uint64>` "Seed Value", with
> `params.linkParameters(k_UseSeed_Key, k_SeedValue_Key, true)`. In
> `executeImpl`, if "Use Seed" is off the seed is taken from
> `std::chrono::steady_clock::now().time_since_epoch().count()`. The seed
> actually used is written into a top-level UInt64 array (created in preflight
> via `CreateArrayAction`, named by "Stored Seed Value Array Name") for
> reproducibility, then passed to `std::mt19937_64`.

> **Units note:** `Physical Size`, `Physical Spacing`, and the `Theta List`
> correlation lengths must share the same length unit. Internally the
> simulation only uses the dimensionless ratio `lag / theta`, so the absolute
> unit is irrelevant *as long as it is consistent*. The MATLAB defaults were
> mm-scale; document this clearly so users do not mix µm size with mm theta.

### 2.2 Outputs

The filter **creates a new** Image Geometry (it does not write into the ODF
geometry):

- **Geometry:** default name **`MTR Microstructure`**; dims
  `n_i = round(Size_i / Spacing_i)`, origin `(0,0,0)`, spacing =
  `Physical Spacing`. Created in preflight via `CreateImageGeometryAction` (all
  inputs are parameters, so dims are known at preflight).
- **Cell arrays:**
  | Array (default name) | Type | Comps | Notes |
  |---|---|---|---|
  | `MTRIds` | Int32 | 1 | Values start at **1** (0 reserved, matches FeatureIds convention). |
  | `Eulers` | Float32 | 3 | Bunge `phi1, PHI, phi2` [radians]. |
  | `Polar Colors` | UInt8 | 3 | RGB. **Created only when** "Generate Polar Coloring" is on. |

Downstream, users can run the stock **Compute IPF Colors** and **Write Image**
filters for additional visualization; only the bespoke MATLAB polar coloring is
built in.

---

## 3. Algorithm Flow & Technical Concerns

`MTRSim::operator()` performs:

1. **Read ODF** — convert the selected Float64 cell arrays + ODF geometry into a
   `std::vector<mtrsim::ODFComponent>`: copy `odfVal`, normalize so each
   component sums to 1, and derive `phi1Bins/PHIBins/phi2Bins` (radians) from
   the geometry dims + degree-spacing.
2. **Build `SimulationParams`** from the filter parameters (Size→`xLen/yLen/zLen`,
   Spacing→`dx/dy/dz`, `volumeFractions`, `thetaList`, `seed`).
3. **PGRF** — `PGRFSimulation::run` → 1-based `mtrIndex` per voxel.
4. **Sample orientations** — `ODFSampler::sampleN` per component against a
   uniform reference ODF.
5. **Assign** — per voxel, pick the orientation from the component named by
   `mtrIndex`.
6. **Write** MTR Index + Euler arrays into the new geometry's cell arrays.
7. **Polar color** (optional) — `IPFMapper::eulerToColors(..., MatLab)` →
   UInt8 RGB.

### Concerns that drive correctness

1. **Voxel index remapping (highest risk).** SIMPLNX requires cell data laid
   out **`z` (slowest) → `y` → `x` (fastest)** in memory — index
   `(z·ny + y)·nx + x`. The standalone driver instead iterates `z → x → y`
   (`main.cpp`), producing `k = ((z)·nx + x)·ny + y` — the column-major MATLAB
   ordering. The algorithm must remap the simulation's `z,x,y` output into the
   SIMPLNX `z,y,x` layout when filling cell arrays, or the field comes out
   transposed. This remap gets a dedicated small-grid deterministic test.
2. **`buildUniformODF()` exposure.** Currently in `main.cpp`'s anonymous
   namespace, hardcoded to 72×36×72. Move it into `LibMTRSim` and
   **parameterize by grid dims** so the uniform reference is derived from the
   actual ODF geometry rather than assuming a 5° grid.
3. **Units consistency** — see the units note in §2.1.
4. **ODF→component helper** — a new small building block (Float64 arrays + ODF
   geometry → `ODFComponent`), the inverse of AJ's `ReadMTRSimODF` write path.

---

## 4. Testing Strategy

Philosophy: `LibMTRSim` already has its own unit tests (`tests/`) proving the
simulation's numerical correctness. The **filter** tests target only what the
filter adds on top of the library.

### 4.1 Component-level, exact / deterministic
- ODF read-back: a known ODF geometry → expected `ODFComponent` (bin centres,
  normalized values).
- `buildUniformODF` bin-centre values for given grid dims.
- Voxel index remap on a tiny grid (e.g. 2×3×2) — exact positional check.
- Polar-color LUT: a known Euler triple → expected RGB.

### 4.2 End-to-end, statistical
- Fixed seed + the **same ODF exemplar and inputs used by the LibMTRSim test**,
  so the filter's results are compared against the *same reference data*. This
  proves the filter wired the library up correctly.
- Assertions use **remap-invariant** statistics: empirical volume fractions ≈
  targets (tolerance), MTR-index value set `{1..N}`, Euler ranges valid,
  per-component mean orientation near the ODF peak. No positional / bit-exact
  array comparison (would be fragile across the 3 CI platforms).

### 4.3 Preflight / error tests
- Volume-fraction column count ≠ `numComponents`.
- Theta rows `< numComponents - 1`.
- Volume fractions do not sum to 1.
- Empty / invalid ODF component selection.

### 4.4 Exemplar data
Stored compressed in-repo, consistent with the AJ convention.

---

## 5. CI (already added — validate)

`.github/workflows/` already contains `linux.yml`, `macos.yml`, `windows.yml`,
`format_pr.yml`, `format_push.yml`, plus issue/PR templates — modeled on
`SimplnxReview`. Each build job clones `simplnx`, configures with
`-DSIMPLNX_EXTRA_PLUGINS="MTRSim" -DSIMPLNX_PLUGIN_ENABLE_MTRSim=ON
-DSIMPLNX_MTRSim_SOURCE_DIR=<workspace>`, builds, and runs `ctest`. vcpkg
binary caching comes from the BlueQuartz NuGet package registry.

Remaining work:
- Confirm the **embedded LibMTRSim vcpkg dependencies** (Eigen, spdlog, CLI11,
  nlohmann-json, hdf5, stb) resolve inside the simplnx build on all three
  platforms.
- Commit + push and confirm the workflows trigger and pass green.

---

## 6. File Plan

New / changed files:

```
src/MTRSim/Filters/MTRSimFilter.{hpp,cpp}          # filter
src/MTRSim/Filters/Algorithms/MTRSim.{hpp,cpp}     # algorithm
src/LibMTRSim/...                                   # move/parameterize buildUniformODF;
                                                    #   add ODF-geometry→ODFComponent helper
src/MTRSim/MTRSimPlugin.cpp                         # register MTRSimFilter
test/MTRSimTest.cpp                                 # filter unit + statistical tests
test/CMakeLists.txt                                 # add test
docs/.../MTRSimFilter.md                            # filter documentation (Milestone AL polishes)
.github/workflows/*.yml                             # already added — validate green
```

---

## 7. Out of Scope (this milestone)

- Self-paced tutorial, example pipelines, and final documentation polish →
  Milestone AL.
- Any change to the MATLAB reference code or the ODF HDF5 on-disk format.
- Exposing `nuggetVariance` (unused by the simulation).

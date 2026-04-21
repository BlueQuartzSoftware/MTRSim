# Milestone AJ Design — MTR Virtual Representations Task 2

**Date:** 2026-04-21
**Milestone:** AJ — Design and Development of Supporting MTR Representation Generation Codes
**Repository:** MTRSim (single-repo approach; existing library repo is promoted to host the new simplnx plugin)

---

## 1. Purpose and Context

Milestone AJ delivers the data-model plumbing needed for MTR virtual-representation workflows inside DREAM3D-NX. The forward simulation algorithms (Milestone AH, delivered) live in `libmtrsim`. AJ adds the SIMPLNX-framework integration: filters that import, export, and construct the Orientation Distribution Functions (ODFs) that feed the Milestone AK simulation filter.

The work is purely additive to simplnx core — no core changes — and preserves byte-exact round-trip compatibility with the existing MATLAB-produced HDF5 ODF format so that ongoing MATLAB-vs-C++ validation continues to work.

### SOW Deliverable (verbatim)

> New or extended data structures and supporting code shall be delivered that enable the storage and manipulation of microstructure descriptors relevant to MTR workflows within the DREAM3D-NX platform.

### SOW Acceptance Criteria (verbatim)

> The enhancements shall be integrated into the DREAM3D-NX framework and verified to support access, modification, and serialization of relevant descriptor data without disrupting existing workflows.

### Scope summary

AJ comprises three concrete SIMPLNX filters, two new library helpers, a plugin + dual-build CMake setup, and a test suite. All are delivered as a single coherent milestone.

---

## 2. Repository and Plugin Layout

The existing `MTRSim` repository is promoted to host both the library and the new simplnx plugin. Single repository; no split.

```
MTRSim/
├── CMakeLists.txt                     ← top-level: builds the simplnx plugin (always)
├── cmake/
├── src/
│   ├── libmtrsim/                     ← existing library (unchanged structure)
│   │   ├── <existing sources>
│   │   ├── ODFBuilder.{hpp,cpp}       ← NEW helper
│   │   └── SymmetricEulers.{hpp,cpp}  ← NEW helper
│   └── MTRSim/                        ← NEW: simplnx plugin source tree
│       ├── MTRSimPlugin.{hpp,cpp}
│       └── Filters/
│           ├── ImportMTRSimODFFilter.{hpp,cpp}
│           ├── ExportMTRSimODFFilter.{hpp,cpp}
│           ├── ComputeODFFromEulerAnglesFilter.{hpp,cpp}
│           └── Algorithms/
│               ├── ImportMTRSimODF.{hpp,cpp}
│               ├── ExportMTRSimODF.{hpp,cpp}
│               └── ComputeODFFromEulerAngles.{hpp,cpp}
├── app/                                ← existing CLI driver (retained; opt-in build)
├── tests/
│   ├── <existing library tests>
│   └── plugin/                         ← NEW: filter tests
├── wrapping/python/                    ← existing pybind11 bindings (retained; opt-in build)
├── configs/, data/, matlab/, tools/    ← unchanged
└── docs/
    └── superpowers/specs/              ← this document
```

The plugin file/directory layout (`src/MTRSim/Filters/{Name}Filter.{hpp,cpp}` and `Filters/Algorithms/{Name}.{hpp,cpp}`) matches what the simplnx `scripts/make_filter.py` generator produces. That generator will be used to scaffold each filter and its test skeleton.

### Build modes

- **Default build:** simplnx `MTRSim` plugin only. `libmtrsim` sources are compiled directly into the plugin target via `target_sources()` — no separately-installed shared library in deployment.
- **Opt-in (`MTRSIM_BUILD_STANDALONE_LIB=ON`):** additionally builds `libmtrsim.{so,dylib,dll}`, the `mtrsim` CLI (`app/`), and Python bindings (`wrapping/python/`). Preserves the existing MATLAB-comparison and scripting workflows.

### Dependencies

- Plugin: `EbsdLib::EbsdLib`, `HDF5::HDF5`, simplnx core libs declared by `make_filter.py`.
- Standalone lib+CLI: adds `spdlog`, `CLI11`, `nlohmann_json`.
- Plugin has no runtime dependency on `libmtrsim`-as-a-shared-library.

CMake file authoring is owned by the project lead; this spec documents the target shape and options.

---

## 3. Extended Data-Structure Usage Convention: "ODF as ImageGeom"

This is the load-bearing design decision that shapes every filter below. Recorded here once, reused throughout.

**An ODF in the SIMPLNX DataStructure is represented as an `ImageGeom` with one or more Float64 single-component cell-data arrays, one array per ODF component.**

| Property | Value |
|---|---|
| Geometry type | `ImageGeom` |
| Default geometry name | `ODF` |
| Dimensions `(X, Y, Z)` | `(num phi2 bins, num PHI bins, num phi1 bins)` |
| Spacing `(X, Y, Z)` | `(phi2 bin size, PHI bin size, phi1 bin size)` — **in degrees** |
| Origin | `(0.0, 0.0, 0.0)` |
| Cell-data arrays | N × `Float64` single-component, one per component |
| Default array names (on import) | `component_0`, `component_1`, … (from source HDF5 group names) |
| Default array names (on construction) | User-supplied, default `"Component 1"` for first-run |

### Axis / index mapping

MATLAB convention (and the existing HDF5 format) uses:
- `phi1` slowest-varying (maps to Z)
- `PHI`  middle      (maps to Y)
- `phi2` fastest-varying (maps to X)

The row-major linear index formula (C order) matches the existing MATLAB / C++ convention:
```
ix = (i_phi1) * (nPHI * nphi2) + (i_PHI) * (nphi2) + (i_phi2)
```

### Known UX consideration

Standard ODF plots slice the phi1–PHI plane at fixed phi2 values — this is the "textbook" orientation. With our axis mapping, that corresponds to DREAM3D-NX's **YZ image-slicing plane** at a chosen X slice. Each filter's user-facing documentation must explain this mapping clearly.

**Future enhancement (out of AJ scope):** an optional "Transpose ODF Axes" utility filter that re-lays the arrays into the textbook orientation for users who prefer it. AJ filters always use the MATLAB-compatible mapping so round-trip is preserved.

### Units

- Degrees in all filter parameters, in the ImageGeom spacing, and in the displayed UI.
- Radians only inside algorithm code (conversion at entry/exit points).
- No radians cross filter or DataStructure boundaries.

### Crystal symmetry

Not stored on the ODF geometry. When needed (EBSD → ODF build in AJ, MTR simulation in AK), the user runs the existing "Create Ensemble Info" filter upstream and passes the resulting **Crystal Structures** (ensemble) and **Phases** (cell) arrays as separate filter parameters.

### Bin-array regeneration

The traditional `phi1_bins` / `PHI_bins` / `phi2_bins` edge arrays are **not stored** in the DataStructure — they're regenerable on-demand from the geometry's spacing and dims. The `ImageGeom` API provides coordinate-to-index and index-to-coordinate helpers directly. Bin arrays are only reconstituted during export, when round-tripping back to the MATLAB HDF5 format that requires them.

---

## 4. Filter: `ImportMTRSimODFFilter`

**Purpose:** Read a MATLAB-compatible MTRSim ODF HDF5 file into an `ImageGeom` + N Float64 cell-data arrays.

### Parameters

| Key | Type | Default | Notes |
|---|---|---|---|
| `input_file` | `FileSystemPathParameter` | — | `.h5` / `.hdf5` only; must exist |
| `output_image_geometry` | `DataGroupCreationParameter` | `/ODF` | Path where the new ImageGeom is created |
| `cell_attribute_matrix_name` | `DataObjectNameParameter` | `Cell Data` | Standard SIMPLNX convention |

Component array names are derived automatically from the HDF5 group names (`component_0`, `component_1`, …). Users may rename them downstream via the standard "Rename Data Object" filter.

### Preflight

1. Open input file with HDF5.
2. Read `/ODF_best/num_components` (int64).
3. For each `component_0 .. component_{N-1}` group, read the sizes of `ODFval`, `phi1_bins`, `PHI_bins`, `phi2_bins`, and load the full bin arrays (cheap: a few hundred doubles).
4. **Strict validation** — any failure is a preflight error, not a warning:
   - Components present contiguously from `0` to `N-1` (no gaps).
   - `ODFval` size equals `(phi1_bins.size - 1) × (PHI_bins.size - 1) × (phi2_bins.size - 1)` for every component.
   - `phi1_bins`, `PHI_bins`, `phi2_bins` are byte-exact identical across all N components. Any mismatch reports the offending component index.
   - Bin arrays are monotonically increasing and uniform-stepped (sanity check).
5. Derive `ImageGeom` spacing from bin-array step (converted rad → deg).
6. Emit `CreateImageGeometryAction` + N `CreateArrayAction` (Float64, single-component, tuple shape matches derived dims).

### Execute

1. Re-open the file.
2. For each component, copy `ODFval` into the corresponding cell-data array (layout already matches row-major flat; no reshape required).
3. Close file.

### Shared library helper

```
namespace mtrsim {

struct ODFFileMetadata {
  int64_t numComponents;
  std::array<int64_t, 3> dimsZYX;             // phi1, PHI, phi2 bin counts
  std::array<double, 3>  spacingDegZYX;       // phi1, PHI, phi2 bin size in degrees
};

struct ODFFileComponent {
  std::vector<double> values;                  // row-major ODFval, size = prod(dims)
};

ODFFileMetadata readODFMetadata(const std::filesystem::path& file);
std::vector<ODFFileComponent> readODFComponents(const std::filesystem::path& file);

}  // namespace mtrsim
```

Preflight uses the metadata-only reader. Execute uses the full reader. Both live in `libmtrsim` and are unit-testable without the simplnx plugin loaded.

### Tests

- **Happy path:** round-trip against a known-good exemplar (gzipped copy of `data/simulation_ODF.h5`). Float64 array equality is bit-exact.
- **Negative paths:** malformed file, missing component, bin-array mismatch between components, non-`.h5` extension, `ODFval` size mismatch.

---

## 5. Filter: `ExportMTRSimODFFilter`

**Purpose:** Write a SIMPLNX ODF (ImageGeom + N Float64 cell-data arrays) back to the MATLAB-compatible HDF5 layout. Round-trips losslessly with `ImportMTRSimODFFilter`.

### Parameters

| Key | Type | Default | Notes |
|---|---|---|---|
| `output_file` | `FileSystemPathParameter` | — | Save-file mode; `.h5` / `.hdf5` extension enforced |
| `input_image_geometry` | `GeometrySelectionParameter` | — | Restricted to `IGeometry::Type::Image` |
| `odf_components` | `MultiArraySelectionParameter` | — | User picks which cell arrays to export. Must be Float64, single-component, on the selected geometry, all sharing the same tuple count |

### Preflight

1. Selected geometry must be an `ImageGeom`.
2. Selected arrays must live on the cell attribute matrix of the selected geometry.
3. At least one array selected.
4. All selected arrays share tuple count equal to `dims.x * dims.y * dims.z`.
5. Tuple count along each axis is ≥ 2 (degenerate geometries rejected).

### Execute

1. Create / overwrite the output HDF5 file.
2. Write `/ODF_best/num_components` = **number of selected arrays**.
3. For each selected array `i` (in user-selection order):
   - Create group `component_i`.
   - Write `ODFval` dataset — direct byte copy of the array buffer (already row-major Float64).
   - Generate `phi1_bins`, `PHI_bins`, `phi2_bins` on the fly from the geometry's spacing and dims (N+1 edges for N bins along that axis), convert deg → rad, write as Float64 datasets.
4. Close file.

### Semantics note (documented; not enforced)

The HDF5 output uses `component_0 .. component_{N-1}` where `i` is the user's **selection index**, not a name derived from the source array's DataPath. If the user selects only `component_1` and `component_3` from an imported geometry, the exported file's `component_0` will hold the input's `component_1` data. The MultiArraySelection order is the source of truth. Documented in filter help text.

### Tests

- **Round-trip:** Import → Export → `h5diff` against original exemplar. Expect byte-exact match.
- **Negative paths:** geometry not `ImageGeom`, selected array on wrong geometry, selected array not Float64, zero arrays selected.

---

## 6. Filter: `ComputeODFFromEulerAnglesFilter`

**Purpose:** Build an ODF from an EBSD scan (per-voxel Euler angles + per-voxel phase labels + per-phase crystal structures), applying crystal-symmetry expansion and optional tri-linear neighbor-bin smoothing. Functionally equivalent to MATLAB `calc_ODF.m` + `symmetric_euler_angles.m`.

This is the only AJ filter where we write genuinely new algorithm code. The other two filters are I/O.

### Parameters

| Key | Type | Default | Notes |
|---|---|---|---|
| `output_mode` | `ChoicesParameter` | `"Create New ODF Geometry"` | Choices: "Create New ODF Geometry" \| "Append to Existing ODF Geometry". Drives linked-parameter visibility |
| `apply_smoothing` | `BoolParameter` | `true` | Tri-linear face/edge/corner smoothing (per `calc_ODF.m`) |
| `bin_size_deg` | `Float32Parameter` | `5.0` | Uniform across phi1 / PHI / phi2. Active only in "Create New" mode. Preflight errors if `360 / bin_size` is not a positive integer (phi1/phi2) or `180 / bin_size` is not a positive integer (PHI) |
| `euler_angles` | `ArraySelectionParameter` | — | `Float32`, 3 components, cell-level data |
| `phases` | `ArraySelectionParameter` | — | `Int32`, 1 component, same cell attribute matrix as `euler_angles` |
| `crystal_structures` | `ArraySelectionParameter` | — | `UInt32`, 1 component, ensemble-level data (produced upstream by "Create Ensemble Info") |
| `mask` | `ArraySelectionParameter` (optional) | — | `Bool`, 1 component, same cell attribute matrix as `euler_angles`. When present, only masked-true voxels contribute |
| **Create New mode parameters:** | | | |
| `output_image_geometry` | `DataGroupCreationParameter` | `/ODF` | |
| `cell_attribute_matrix_name` | `DataObjectNameParameter` | `Cell Data` | |
| `component_name` | `DataObjectNameParameter` | `Component 1` | |
| **Append mode parameters:** | | | |
| `existing_odf_geometry` | `GeometrySelectionParameter` | — | Restricted to `IGeometry::Type::Image` |
| `component_name` | `DataObjectNameParameter` | `Component 1` | Must not collide with an existing cell-data array name on the target geometry |

### Algorithm

For each voxel `i` where the mask (if present) is true and the phase is non-zero:

1. Read `(phi1, PHI, phi2)` from `euler_angles[i]`.
2. Look up `crystal_structures[phases[i]]` to get the crystal-symmetry enum.
3. Use `SymmetricEulers` helper to expand to the full list of symmetric-equivalent tuples for that crystal system (via EbsdLib).
4. For each symmetric tuple:
   - Compute bin index `(i_phi1, i_PHI, i_phi2)` from `bin_size_deg`.
   - If `apply_smoothing` is **on**, distribute 1.0 of contribution across:
     - Center bin ×  `0.332`
     - 6 face-neighbors × `0.448 / 6` each
     - 12 edge-neighbors × `0.16 / 12` each
     - 8 corner-neighbors × `0.06 / 8` each
     - All neighbor accumulations per symmetric tuple. Boundary bins wrap via Bunge-angle periodicity (phi1 and phi2 modulo 2π; PHI reflects at 0 and π under crystal symmetry) exactly as implemented in `matlab/calc_ODF.m` (see the `jf_minus` / `jf_plus` / `kf_minus` / `kf_plus` / `lf_minus` / `lf_plus` index logic).
   - If `apply_smoothing` is **off**, add `1.0` to the center bin only.
5. Normalize the entire `ODFval` array by `N` = count of voxels contributing (not the expanded count). Matches MATLAB normalization.

### New library helpers (in `libmtrsim`)

Both helpers live in `src/libmtrsim/`, are namespaced under `mtrsim::`, and ship with library-level unit tests independent of the plugin.

```
// SymmetricEulers.hpp
namespace mtrsim {
  // Returns the list of symmetric-equivalent Euler tuples (Bunge phi1/PHI/phi2, radians)
  // for the given input and crystal system. Thin wrapper over EbsdLib's orientation
  // operator classes. Result length depends on the crystal system (e.g., HCP = 12).
  std::vector<std::array<double, 3>>
  expandSymmetric(double phi1, double PHI, double phi2, uint32_t ebsdLibCrystalCode);
}

// ODFBuilder.hpp
namespace mtrsim {
  struct ODFBuildParams {
    int32_t nphi1;                 // bins along phi1 (ImageGeom Z)
    int32_t nPHI;                  // bins along PHI  (ImageGeom Y)
    int32_t nphi2;                 // bins along phi2 (ImageGeom X)
    double  binSizeDeg;            // uniform across all three axes
    bool    smoothing;             // tri-linear face/edge/corner distribution
  };

  // Accumulates Euler-tuple contributions into an ODFval vector, applying the
  // smoothing distribution if enabled. Input tuples are in radians.
  // `values` is pre-allocated to nphi1*nPHI*nphi2, initially zeros.
  void accumulate(const std::vector<std::array<double, 3>>& symmetricEulers,
                  const ODFBuildParams& params,
                  std::vector<double>& values);

  // Final normalization pass: divides `values` by `normalizer`.
  void normalize(std::vector<double>& values, double normalizer);
}
```

**EbsdLib policy:** All orientation math (symmetry operators, Euler↔matrix conversions, crystal-system enumeration) routes through EbsdLib. The spec and future PRs will document any divergence if genuinely unavoidable.

### Preflight

1. Mode-gated parameter visibility validation.
2. `euler_angles`, `phases`, `mask` (if present) all live on the same cell attribute matrix.
3. `crystal_structures` lives on an ensemble attribute matrix.
4. `bin_size_deg` produces integer bin counts for both 360° axes and the 180° axis.
5. **Append mode only:**
   - `existing_odf_geometry` spacing is uniform across all three axes (non-uniform bin sizes are unsupported).
   - `existing_odf_geometry` spacing matches the implied bin size, OR the spec uses the existing spacing as source-of-truth and `bin_size_deg` is hidden/ignored. **Decision: existing geometry's spacing wins; `bin_size_deg` is grayed out in Append mode.**
   - `component_name` does not already exist as a cell-data array on the target geometry.
6. Emit `CreateImageGeometryAction` + `CreateArrayAction` (Create New mode), or just `CreateArrayAction` (Append mode).

### Parallelization

The output array is small (186,624 × 8 bytes ≈ 1.5 MB at 5° resolution). Use the standard simplnx `ParallelDataAlgorithm` pattern with per-thread local accumulators merged at the end. No atomics needed; no contention. Follow the `bluequartz-skills:thread-safety` and `bluequartz-skills:progress-messaging` conventions.

### Tests (tolerance-based)

- **Ground truth:** a small EBSD fixture (≤1000 points, HCP phase) run through MATLAB `calc_ODF.m` with known fixed smoothing settings; resulting `ODFval` saved as a compressed HDF5 reference.
- **Tolerance:** target `max(|diff|) < 1e-10`. Actual threshold tuned when tests are implemented, depending on observed FP summation-order differences.
- **Test matrix:**
  - Smoothing on / off.
  - Single-phase vs two-phase dataset with mask selecting one phase.
  - Create New mode vs Append mode.
  - HCP symmetry (primary use case); cubic as a sanity check if fixture available.
- **Negative paths:** non-integer bin-size divisor, component-name collision in Append mode, missing ensemble data, `phases` and `euler_angles` on different cell attribute matrices.

---

## 7. Testing Infrastructure

### Exemplar data storage location — **DEFERRED**

Options to be decided before the first test-data commit:

| Option | Pros | Cons |
|---|---|---|
| In-repo under `tests/data/` (gzipped HDF5) | Simple; git-tracked; no external dependency | Bloats repo for large EBSD fixtures |
| GitHub release artifacts | Scales to large binaries | Release-management overhead |
| Local webserver (shared simplnx infra) | Matches existing simplnx plugin pattern | Shared URL cross-plugin coupling |

**Recommendation to revisit at testing phase:** In-repo for ODF-only fixtures (<10 MB each). External store for EBSD inputs if needed for the EBSD→ODF tests.

### Exemplar format

Gzipped HDF5 files (HDF5 has built-in gzip dataset compression; no tar layer needed). Re-evaluated alongside the storage-location decision.

### Test tooling

- Catch2 (existing plugin convention).
- `make_filter.py` scaffolds each filter's test skeleton.
- Error-path coverage for every preflight failure mode listed above.

---

## 8. Acceptance-Criteria Mapping

| SOW Criterion | AJ Evidence |
|---|---|
| Integrated into DREAM3D-NX | Plugin builds against simplnx core; filters listed under the "MTRSim" plugin group in DREAM3D-NX UI |
| Access | `ImportMTRSimODFFilter` reads existing MATLAB-format HDF5; geometry + arrays become visible to all downstream filters |
| Modification | `ComputeODFFromEulerAnglesFilter` creates or appends ODF components in-place in the DataStructure |
| Serialization | `ExportMTRSimODFFilter` round-trips to MATLAB-compatible HDF5; byte-exact against imported exemplar |
| Without disrupting existing workflows | Purely additive plugin — no changes to simplnx core; MATLAB workflow preserved via round-trip guarantee |

### Report framing: deliverable adequacy

Although AJ reuses existing SIMPLNX types rather than introducing a bespoke `ODFData` class, the deliverable is satisfied by three concrete categories of artifact:

1. **Extended Data-Structure Usage Convention** — the "ODF as ImageGeom" convention documented in Section 3 is a reusable contract future filters and workflows will consume.
2. **New library helpers** — `ODFBuilder` and `SymmetricEulers` in `libmtrsim` (the "supporting code" per the deliverable language).
3. **Three SIMPLNX filters** — access, modification, and serialization surfaces for MTR descriptor data.

The report narrative will frame this as *"we extended the DataStructure's usage conventions to represent Euler-space distributions without requiring a bespoke type — this is the SIMPLNX-idiomatic approach."* Similar precedents exist in the simplnx codebase where established types model domain-specific concepts via convention rather than new classes.

---

## 9. Out-of-Scope / Future Work

- **Transpose ODF Axes helper filter** (textbook phi1-X / PHI-Y / phi2-Z orientation). Nice-to-have but breaks MATLAB round-trip by default; defer until a concrete user request arises.
- **Non-uniform bin-size support** in `ComputeODFFromEulerAnglesFilter`. Current design enforces uniform bins across all three axes.
- **Additional EBSD reader variants** beyond what EbsdLib and existing simplnx filters already provide. AJ's EBSD→ODF filter consumes Euler-angle cell arrays produced by the existing simplnx EBSD-read filters; it does not re-invent file readers.
- **MTR Index renumbering to zero-based.** MATLAB / C++ / AK retain 1-based indexing (where 0 is reserved for "unindexed" per DREAM3D convention).
- **AK / AL milestones.** Covered in separate design specs.

---

## 10. Deferred Decisions

Explicit TBDs to revisit before or during implementation:

| Item | When to decide |
|---|---|
| Exemplar storage location (in-repo vs GitHub releases vs local webserver) | First test-data commit |
| Exact tolerance threshold for `ComputeODFFromEulerAnglesFilter` tests | When first reference comparison is run |
| Whether the Append-mode filter's `bin_size_deg` parameter is grayed out (UI convention) or hidden (structural) | During filter-parameter UI construction |

---

## 11. High-Level Implementation Order

Informational sequencing for the writing-plans step; not binding. All of the below are within the single AJ milestone and can be PR'd separately.

1. Plugin scaffolding + dual-build CMake (`MTRSIM_BUILD_STANDALONE_LIB` option, plugin target pulls in `libmtrsim` sources).
2. `ImportMTRSimODFFilter` + `readODFMetadata` / `readODFComponents` library helpers + tests.
3. `ExportMTRSimODFFilter` + round-trip tests against #2's exemplar.
4. `SymmetricEulers` + `ODFBuilder` library helpers + library-level unit tests.
5. `ComputeODFFromEulerAnglesFilter` + tolerance-based tests using MATLAB reference.
6. Report draft, acceptance-criteria mapping, deliverable framing.

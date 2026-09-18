# Combined ODF Computation and Section Presentation Design

## Status

Approved design, revised 2026-09-18. This replaces the earlier design for a
separate `WriteODFSectionsFilter`.

## Context

`ComputeODFFilter` produces a quantitative orientation distribution function
(ODF) as a scalar array on a regular Bunge-Euler Image Geometry. The initial
implementation added a second public filter to convert that grid into a
conventional constant-phi2 section montage.

That separation created two avoidable user hazards:

1. The ODF units were only a Compute-filter argument and were not retained in
   the Data Structure. The presentation filter therefore asked users to
   remember whether the array contained MUD or Count-Density.
2. The ODF could contain multiple phases while a standard Euler-section
   presentation requires one well-defined Laue class.

The filters have not been released, so the final public design combines
computation and presentation in `ComputeODFFilter`, enforces one selected phase,
and persists unit provenance with the quantitative ODF.

The existing pole-figure implementation remains the architectural precedent:
EbsdLib owns reusable scientific rendering, while MTRSim owns SIMPLNX Data
Structure integration, pipeline actions, and file output.

## Goals

- Make `ComputeODFFilter` the only public ODF computation/presentation filter.
- Compute one phase per execution using an explicit phase index.
- Preserve the quantitative ODF Image Geometry and scalar component array.
- Persist each component's units in serializable Data Structure objects.
- Optionally generate a conventional constant-phi2 montage as PNG, RGB Image
  Geometry, or both.
- Always render the standard presentation in MUD while identifying whether the
  source array was already MUD or converted from Count-Density.
- Use conventional Greek Euler symbols, degree symbols, and readable angular
  tick marks.
- Support pure pipeline preflight where upstream ensemble values are not yet
  materialized.
- Retain the reviewed EbsdLib scientific preparation and compositor APIs.

## Non-goals

- A standalone public `WriteODFSectionsFilter`.
- Backward compatibility for the unreleased standalone filter or its UUID.
- Manual unit selection or override in the presentation workflow.
- Three-dimensional Euler-space rendering.
- Sigma, omega, axis-angle, or arbitrary section families.
- Arbitrary user-supplied phi2 section angles.
- Rendering multiphase ODF values through one Laue class.
- Supporting cropped, transposed, or nonuniform input ODF grids.

## Approved user-facing decisions

| Decision | Choice |
|---|---|
| Public filter | `ComputeODFFilter` only |
| Computation scope | One required phase index per execution |
| Numeric output | Full ODF Image Geometry and scalar component array |
| Numeric output units | User selects MUD or Count-Density |
| Persistent units | One-value StringArray per ODF component |
| Presentation units | MUD |
| Count-Density presentation label | `MUD (from Count-Density)` |
| Default section count | 6 |
| Section placement | Evenly spaced over the Laue-class phi2 range; upper endpoint excluded |
| Default layout | 3 columns by 2 rows |
| Default color table | `Rainbow Blended White` |
| Scale control | Automatic or Manual MUD range |
| Presentation outputs | Optional PNG and optional RGB Image Geometry |
| Section title | Greek phi2 with subscript and degree symbol, for example `φ₂ = 10°` |
| Axis titles | `φ₁` and `Φ` |
| Axis ticks | 10 degrees when spacing permits, otherwise 20 degrees; numeric labels have no unit suffix |

## Architecture

### Responsibility boundary

EbsdLib owns:

1. Full-grid validation.
2. Count-Density-to-MUD conversion.
3. Laue-class Euler plotting limits.
4. Exact periodic phi2 section interpolation.
5. Common scaling and color control interpolation.
6. Panel layout and annotated RGBA composition.
7. Greek section/axis typography, ticks, labels, and shared color bar.

MTRSim owns:

1. SIMPLNX filter parameters and parameter linking.
2. Phase and optional-mask selection.
3. Numeric ODF accumulation and output actions.
4. Persistent per-component units metadata.
5. Adapting SIMPLNX arrays to EbsdLib.
6. Optional RGB Image Geometry population.
7. Atomic PNG output.

The public filter is combined. Internally, presentation remains a focused
helper/algorithm so ODF accumulation is not mixed with canvas composition in
one large source file.

### EbsdLib components

The reviewed EbsdLib components remain:

- `ODFSectionUtilities`
  - Grid validation.
  - Count-Density-to-MUD conversion.
  - Laue plotting limits.
  - Exact section-angle generation.
  - Periodic cell-center interpolation along phi2.
- `ODFSectionChrome`
  - Section titles, axes, ticks, shared color bar, title, and phase information.
- `ODFSectionCompositor`
  - Scientific preparation orchestration, colorization, layout, and RGBA output.

No SIMPLNX type, path, parameter, or preset name enters EbsdLib.

### MTRSim internal presentation helper

The existing standalone-filter algorithm is retained only as an internal
presentation helper after being renamed or refactored to remove any dependency
on `WriteODFSectionsFilter`.

It consumes:

- The computed ODF array and Image Geometry.
- The selected phase's crystal structure and material name.
- The known numeric output units from the same filter execution.
- Presentation parameters and output paths.

It produces optional PNG/RGB outputs and returns contextual errors to
`ComputeODFFilter`.

## Single-phase computation contract

`ComputeODFFilter` requires:

- Cell Euler Angles.
- Cell Phases.
- Crystal Structures.
- Material Names.
- Phase Index, default 1.
- Optional Bool/UInt8 mask.

A voxel contributes only when:

```text
phases[voxel] == selectedPhase
AND (mask is disabled OR mask[voxel] is true)
```

Phase zero and every non-selected phase are skipped. The selected phase index
must exist in both ensemble arrays and must contain a supported EbsdLib crystal
structure. Material Names and Crystal Structures must have equal tuple counts.

For multiphase EBSD data, users run the filter once per phase. Append mode can
place multiple phase-specific ODF components on one compatible Euler grid.

This contract guarantees that the numeric ODF, Laue plotting region, material
label, and presentation image describe the same phase.

## Persistent units contract

### Data Structure layout

For every ODF component array, `ComputeODFFilter` creates:

```text
<ODF Image Geometry>/ODF Metadata/<Component Name> Units
```

`ODF Metadata` is a DataGroup. `<Component Name> Units` is a one-tuple
StringArray whose only accepted values are:

- `MUD`
- `Count-Density`

The metadata object is separate from the cell Attribute Matrix so its one-tuple
shape does not conflict with the Euler-grid tuple shape. StringArray/DataGroup
objects are serialized by DREAM3D file writing, unlike in-memory
`DataObject::Metadata`.

### Preflight behavior

The selected output-unit argument is known during preflight. A
`CreateStringArrayAction` therefore initializes the units value in both
preflight and execution modes. Downstream filters in the same pipeline can read
the provenance without waiting for ODF execution.

Create-New mode creates the metadata group and unit array with the geometry.
Append mode reuses or creates `ODF Metadata` and adds the selected component's
unit array. A conflicting metadata object or existing component remains an
error.

No manual units override exists. Missing or invalid units metadata is a
contract error for any future consumer.

## Scientific conventions

### Full-Bunge grid

The numeric ODF covers:

- phi1: `[0 degrees, 360 degrees)`
- PHI: `[0 degrees, 180 degrees)`
- phi2: `[0 degrees, 360 degrees)`

Image Geometry mapping:

| Image axis | Euler axis | Array-order role |
|---|---|---|
| X | phi2 | Fastest varying |
| Y | PHI | Middle |
| Z | phi1 | Slowest varying |

Full-source flat indexing is:

```text
linear = phi1 * (nPHI * nphi2) + PHI * nphi2 + phi2
```

Grid validity uses the same quotient-based Float32-safe rule in computation and
presentation:

```text
abs(360 / spacing - nphi1) <= 1e-6
abs(180 / spacing - nPHI)  <= 1e-6
abs(360 / spacing - nphi2) <= 1e-6
```

### Count-Density-to-MUD conversion

The numeric output preserves the user's selected units. Presentation always
uses MUD. Count-Density is converted before section extraction:

```text
MUD(i, j, k) = CountDensity(i, j, k)
             * 8 * pi^2
             / (dphi1 * dPHI * dphi2 * sin(PHI_center(j)))
```

Angular spacings are radians in this formula. MUD output bypasses conversion.

The shared color-bar title is:

- `MUD` when the numeric output is MUD.
- `MUD (from Count-Density)` when conversion was applied.

The information block also records the numeric output units and rendered units.

### Rectangular Euler plotting regions

phi1 remains `[0 degrees, 360 degrees)` for all Laue classes:

| EbsdLib crystal structure | PHI maximum | phi2 maximum |
|---|---:|---:|
| Hexagonal High (`6/mmm`) | 90 degrees | 60 degrees |
| Cubic High (`m-3m`) | 90 degrees | 90 degrees |
| Hexagonal Low (`6/m`) | 180 degrees | 60 degrees |
| Cubic Low (`m-3`) | 90 degrees | 180 degrees |
| Triclinic (`-1`) | 180 degrees | 360 degrees |
| Monoclinic (`2/m`, EbsdLib unique-axis convention) | 90 degrees | 360 degrees |
| Orthorhombic (`mmm`) | 90 degrees | 180 degrees |
| Tetragonal Low (`4/m`) | 180 degrees | 90 degrees |
| Tetragonal High (`4/mmm`) | 90 degrees | 90 degrees |
| Trigonal Low (`-3`) | 180 degrees | 120 degrees |
| Trigonal High (`-3m`) | 90 degrees | 120 degrees |

Cubic classes use the conventional rectangular Euler plotting box rather than
the nonrectangular strict fundamental zone.

### Section angles

For `N` sections and symmetry maximum `phi2Max`:

```text
phi2Section(s) = s * phi2Max / N,  s = 0 ... N - 1
```

The upper endpoint is excluded. Six Hexagonal High sections are therefore 0,
10, 20, 30, 40, and 50 degrees. Each exact section is periodically
interpolated between neighboring cell-centered phi2 planes.

### Typography and ticks

Section titles use embedded-font Unicode text:

```text
φ₂ = 10°
```

Axis titles are `φ₁` horizontally and `Φ` vertically.

For each axis, use 10-degree tick intervals when the pixel distance between
ticks is at least 24 pixels. Otherwise use 20-degree intervals. Tick values are
plain numbers such as `0`, `20`, and `40`; they do not repeat `deg` or `°`.
The axis titles and section title establish the angular context.

Ticks include both zero and the rectangular plotting maximum. PHI increases
down the page to match prepared section row order.

## Combined `ComputeODFFilter` surface

### Computation inputs

- Apply Smoothing.
- Bin Size in degrees.
- Output Units: Count-Density or MUD.
- Cell Euler Angles.
- Cell Phases.
- Crystal Structures.
- Material Names.
- Phase Index.
- Optional Mask.

### Numeric output

- Create-New or Append mode.
- ODF Image Geometry.
- Cell Attribute Matrix name.
- Component name.
- Per-component persisted units metadata.

### Optional presentation

Presentation runs when either output is enabled:

- `Write PNG`, default false.
- `Save Presentation as Image Geometry`, default true.

Both may be false when only the quantitative ODF is wanted.

Presentation parameters:

- Figure Title.
- Number of phi2 Sections, default 6.
- Sections Per Row, default 3.
- Section Width, default 512 pixels.
- Color Table Preset, default `Rainbow Blended White`.
- Scale Mode: Automatic or Manual.
- Manual Minimum/Maximum MUD.
- Output PNG path when enabled.
- Output Presentation Image Geometry path when enabled.

The presentation Image Geometry contains a UInt8 three-component RGB array
named `ODF Sections RGB` under `Cell Data`.

## Preflight

Preflight always validates:

- Matching cell tuple counts for Euler Angles and Phases.
- Valid selected phase index and equal ensemble tuple counts.
- Valid bin size/full-grid dimensions.
- Create-New/Append compatibility and component collisions.
- Metadata group/units-array collisions.
- Mask type and tuple count when enabled.
- Presentation count/layout/width/color/scale/path parameters.

When upstream Crystal Structures and Material Names values are materialized,
preflight also validates the crystal code/material and reports exact section
angles, Euler limits, and presentation dimensions.

When upstream phase metadata values are unavailable during pure pipeline
preflight:

- Do not index EmptyDataStore values.
- Use the conservative Triclinic layout envelope for output actions.
- Omit phase-dependent updated values entirely. Do not emit values containing
  `Deferred until execution` because SIMPLNX will not update them later.
- Report only values already known, such as phase index, numeric output units,
  grid dimensions/spacing, requested layout, scale mode, and conservative
  presentation dimensions.

Execution validates the materialized phase metadata and resizes a conservative
presentation geometry to the actual Laue-specific dimensions.

## Execution flow

1. Resolve the selected phase, material name, crystal structure, and optional
   mask.
2. Accumulate only selected-phase voxels.
3. Normalize the ODF and produce the selected numeric units.
4. Write the scalar ODF array.
5. Persist the selected units StringArray.
6. If presentation is disabled, finish.
7. Resolve and normalize the selected color-table controls.
8. Invoke EbsdLib with the known numeric units and selected phase metadata.
9. Resize/copy the optional RGB Image Geometry.
10. Write the optional PNG through `AtomicFile`.

Cancellation is checked between phases. A failed render or PNG encode does not
replace an existing destination.

## Removal of the standalone filter

Remove:

- `WriteODFSectionsFilter.hpp/.cpp`.
- Public filter registration and UUID.
- Standalone filter documentation.
- Standalone filter-specific test registration.

Refactor useful execution code into an internal presentation helper consumed by
`ComputeODFFilter`. Tests may remain in a focused presentation test file, but
they exercise the single public Compute filter.

## Error behavior

- Selected phase index zero/out-of-range is an error.
- Unknown crystal structure is an error for the selected phase.
- No contributing selected-phase voxels is an error when presentation is
  requested; numeric-only mode retains the existing all-zero behavior unless
  otherwise validated by the computation contract.
- Negative or non-finite ODF values are errors before rendering.
- All-zero automatic presentation scale is an error.
- A grid with no cell centers inside its Laue plotting limits is an error.
- Missing, conflicting, or invalid units metadata is an internal contract error.
- Output errors include complete paths and underlying writer diagnostics.

## Testing

### EbsdLib

- Greek section-title and axis-title text helpers produce `φ₂`, `φ₁`, `Φ`, and
  `°` exactly.
- Tick interval selection chooses 10 degrees at or above 24 pixels per tick and
  20 degrees below that threshold.
- Tick values contain no unit suffix.
- MUD and Count-Density render with `MUD` and
  `MUD (from Count-Density)` color-bar titles, respectively.
- Existing scientific preparation, axis-order, fractional-spacing, empty-crop,
  scale, color, layout, and pixel tests remain green.

### MTRSim

- Phase selection excludes all non-selected phases without requiring a phase
  mask.
- Optional mask further restricts the selected phase.
- Create-New and Append modes create correct per-component units StringArrays.
- Units metadata values are readable during preflight and survive
  WriteDREAM3D/ReadDREAM3D round-trip.
- Presentation uses the same execution's numeric-unit choice with no manual
  units parameter.
- Default color preset is `Rainbow Blended White`.
- PNG-only, RGB-only, both, and numeric-only modes.
- Conservative pure-pipeline preflight omits all phase-dependent updated values
  and contains no `Deferred until execution` text.
- Materialized preflight reports exact phase/Laue/angle values.
- Count-Density presentation is labeled `MUD (from Count-Density)`.
- Existing atomicity, cancellation, Float32/Float64, fractional-grid, and
  pixel-parity tests remain green.
- Filter list contains no standalone `WriteODFSectionsFilter`.

### Real-data validation

Use `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1.ang`:

- Phase 1 contributes 2,493,822 indexed points.
- Six Hexagonal High sections are 0 through 50 degrees.
- PNG and RGB Image Geometry are byte-identical in RGB channels.
- ODF/presentation values are finite.
- Applied MUD maximum matches independent section extraction.
- Panel feature correlations against MTEX remain in the established
  0.9730–0.9783 range.
- Regenerate the documentation montage using `Rainbow Blended White` and the
  corrected typography/ticks.

## Documentation

- Expand `ComputeODFFilter.md` with phase selection, units metadata,
  presentation parameters, images, and unit-label semantics.
- Remove `WriteODFSectionsFilter.md`.
- Replace the canonical montage with the final combined-filter output.
- Update example/validation pipelines to contain the combined filter only.

## Cross-repository delivery

1. Update EbsdLib chrome/tests for Greek typography, adaptive ticks, and
   unit-aware color-bar titles.
2. Refactor MTRSim to the single public filter and add persisted units metadata.
3. Build with `NX-Com-Qt69-Vtk96-Rel-EbsdLib` and run all EbsdLib/MTRSim tests.
4. Rerun the external `1.ang` validation and update documentation artifacts.
5. Publish an EbsdLib release containing the final renderer API before testing
   the normal `NX-Com-Qt69-Vtk96-Rel` dependency path.

No simplnx OrientationAnalysis source changes are required.

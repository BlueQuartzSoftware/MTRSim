# Standard ODF Section Presentation Design

## Status

Approved design for a reusable EbsdLib ODF-section renderer and an MTRSim
SIMPLNX filter that exposes it.

## Context

`ComputeODFFilter` produces a quantitative orientation distribution function
(ODF) as a scalar array on a regular Bunge-Euler Image Geometry. Although this
representation is appropriate for storage and downstream computation, it does
not directly produce the conventional constant-phi2 section montage used by
microscopists and texture analysts.

The existing pole-figure implementation provides the architectural precedent:
EbsdLib owns the scientific rendering and composite-image generation, while a
SIMPLNX filter owns Data Structure validation, parameter handling, file output,
and optional Image Geometry output.

## Goals

- Add a reusable EbsdLib API for producing a conventional phi2-section ODF
  montage from a regular full-Bunge grid.
- Display phi1 horizontally, PHI vertically, and exact constant-phi2 sections.
- Use the selected phase's Laue class to choose the rectangular Euler plotting
  region and phi2 section range.
- Support MUD and Count-Density inputs while always presenting the result in
  MUD.
- Allow users to select any registered SIMPLNX RGB color-table preset.
- Produce a PNG file, an RGB SIMPLNX Image Geometry, or both.
- Provide deterministic, independently testable scientific transformations and
  rendering.

## Non-goals

- Three-dimensional Euler-space rendering.
- Sigma, omega, axis-angle, or arbitrary section families.
- User-supplied lists of arbitrary phi2 section angles.
- ODF calculation from Euler angles; `ComputeODFFilter` remains responsible for
  that operation.
- Inferring which phase generated an ODF.
- Supporting cropped or nonuniform Euler grids.
- Adding the SIMPLNX wrapper to OrientationAnalysis in this version.

## Approved user-facing decisions

| Decision | Choice |
|---|---|
| Canonical representation | Constant-phi2 Euler sections |
| Default section count | 6 |
| Section placement | Evenly spaced over the Laue-class phi2 range; upper endpoint excluded |
| Default layout | 3 columns by 2 rows |
| Layout configurability | User-selectable sections per row |
| Input units | MUD or Count-Density |
| Rendered units | MUD |
| Crystal symmetry source | Required `Crystal Structures` ensemble array plus phase index |
| Material label source | Required `Material Names` ensemble array |
| Color selection | SIMPLNX `CreateColorMapParameter`, default `Rainbow` |
| Automatic scale | 0 MUD through the maximum displayed value |
| Manual scale | User-supplied minimum and maximum MUD |
| Output | Optional PNG and optional RGB Image Geometry; at least one required |
| SIMPLNX plugin | MTRSim |

## Architecture

### Responsibility boundary

EbsdLib owns all scientific presentation semantics and RGBA composition:

1. Validate the regular Euler grid.
2. Convert Count-Density to MUD when requested.
3. Determine the rectangular Euler plotting region for the Laue class.
4. Generate exact section angles.
5. Interpolate the cell-centered grid to those section angles.
6. Apply the common color scale and color table.
7. Draw the section panels, axes, labels, title, information block, and shared
   color bar.

MTRSim owns all SIMPLNX integration:

1. Resolve Data Paths and validate Data Structure relationships.
2. Read the selected crystal structure and material name.
3. Resolve the SIMPLNX color-table preset.
4. Adapt Float32 or Float64 input into the EbsdLib data contract.
5. Call EbsdLib and translate its diagnostics into SIMPLNX results.
6. Populate the optional RGB Image Geometry.
7. Write the optional PNG atomically.

Neither layer recomputes an ODF from Euler angles.

### EbsdLib components

Add these public utilities under `Source/EbsdLib/Utilities`:

- `ODFSectionUtilities`
  - Grid validation.
  - Count-Density-to-MUD conversion.
  - Laue-class Euler plotting limits.
  - Section-angle generation.
  - Periodic cell-center interpolation along phi2.
- `ODFSectionChrome`
  - Axes, tick labels, section titles, shared MUD color bar, figure title, and
    phase-information block.
- `ODFSectionCompositor`
  - Orchestration, panel colorization, layout, and final RGBA extraction through
    `canvas_ity`.

This separation follows `PoleFigureCompositor` and `PoleFigureChrome` while
keeping the scientific transformations independently testable without image
composition.

### EbsdLib configuration contract

The public configuration contains:

- A non-owning `DoubleArrayType*` containing the ODF scalar values.
- Explicit `(nphi1, nPHI, nphi2)` dimensions.
- Euler-grid origin and spacing in degrees.
- Input-unit enumeration: `MUD` or `CountDensity`.
- EbsdLib crystal-structure/Laue index.
- Section count.
- Sections per row.
- Section width in pixels.
- Automatic or manual scale mode and manual limits.
- Normalized RGB color control points supplied by the caller as repeated
  `[position, red, green, blue]` values.
- Figure title, phase number, material name, and informational labels.

The result contains:

- RGBA image data.
- Composite width and height.
- Exact section angles used.
- Applied MUD minimum and maximum.
- Structured error code and message when generation fails.

The API contains no SIMPLNX Data Paths, Image Geometry types, parameter types,
or preset names.

## Scientific conventions

### Full-Bunge grid

The input must cover the full Bunge cube:

- phi1: `[0 degrees, 360 degrees)`
- PHI: `[0 degrees, 180 degrees)`
- phi2: `[0 degrees, 360 degrees)`

The MTRSim/SIMPLNX Image Geometry mapping is:

| Image axis | Euler axis | Array-order role |
|---|---|---|
| X | phi2 | Fastest varying |
| Y | PHI | Middle |
| Z | phi1 | Slowest varying |

All three axes must have the same positive spacing. Dimensions must equal
`360 / spacing`, `180 / spacing`, and `360 / spacing` for phi1, PHI, and phi2,
respectively. The origin must be zero and the quotients must be positive
integers within floating-point tolerance.

Cell values are located at bin centers. For a five-degree grid, the centers are
2.5, 7.5, and so on.

### Count-Density-to-MUD conversion

Count-Density is converted before section extraction:

```text
MUD(i, j, k) = CountDensity(i, j, k)
             * 8 * pi^2
             / (dphi1 * dPHI * dphi2 * sin(PHI_center(j)))
```

Angular spacings are expressed in radians in this formula. Because the full
grid is cell-centered, `PHI_center` never lies exactly on a Bunge pole. MUD
input bypasses this conversion.

The conversion matches `ComputeODFFilter`, so mathematically equivalent MUD
and Count-Density inputs render identically.

### Rectangular Euler plotting regions

phi1 remains `[0 degrees, 360 degrees)` for all Laue classes. EbsdLib owns the
following PHI and phi2 maxima:

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

For cubic classes this is the conventional rectangular Euler plotting box, not
the nonrectangular strict fundamental zone.

### Section angles and interpolation

For `N` requested sections and symmetry maximum `phi2Max`, section `s` uses:

```text
phi2Section(s) = s * phi2Max / N,  s = 0 ... N - 1
```

The redundant `phi2Max` endpoint is excluded. This reproduces the MTEX default;
for six hexagonal sections the angles are 0, 10, 20, 30, 40, and 50 degrees.

Each exact section is linearly interpolated between its two neighboring
cell-centered phi2 planes. Interpolation is periodic across 360 degrees, so the
zero-degree section interpolates between the final and first planes. No
interpolation is performed along phi1 or PHI; those cell-centered values are
colorized and scaled into the displayed panel.

### Panel geometry and layout

The user controls section width. Section height is derived from the displayed
angular extents so phi1 and PHI use equal angular scale:

```text
sectionHeight = sectionWidth * PHIMax / 360 degrees
```

The default is six sections with three sections per row. Other section counts
use the requested `Sections Per Row`; the final row may be partially filled.
Unused cells are left blank. A single shared color bar is placed to the right of
the montage.

### Color mapping and scaling

The MTRSim filter uses `CreateColorMapParameter` with `Rainbow` as the default.
It resolves the selected preset through `ColorTableUtilities` and passes plain
normalized RGB control points to EbsdLib.

Automatic scaling uses `[0, maximumDisplayedMUD]`, where the maximum is measured
after conversion and exact section interpolation. An all-zero result is an
error because it has no meaningful automatic ODF scale.

Manual scaling requires finite values with `minimum < maximum`. Values below or
above the range use the first or last color, respectively. The shared color bar
always reports the applied MUD range.

## MTRSim filter

### Name and ownership

- C++ name: `WriteODFSectionsFilter`
- Human name: `Write ODF Section Image`
- Plugin: MTRSim
- Algorithm class: `WriteODFSections`

### Parameters

#### Input ODF

- `ODF Image Geometry`: required Image Geometry.
- `ODF Array`: required scalar Float32 or Float64 cell array on that geometry.
- `Input ODF Units`: `MUD` or `Count-Density`.

#### Phase metadata

- `Crystal Structures`: required UInt32 ensemble array.
- `Material Names`: required StringArray ensemble array.
- `Phase Index`: required positive integer, default 1.

The selected ODF must represent the selected phase alone. A pipeline processing
multiphase EBSD data normally creates that ODF with a phase mask. The filter can
validate the selected metadata tuple but cannot prove ODF phase provenance.
When an imported ODF has no ensemble metadata, users can create the required
arrays with `Create Ensemble Info`.

#### Presentation

- `Figure Title`: string.
- `Number of phi2 Sections`: integer, default 6.
- `Sections Per Row`: integer, default 3.
- `Section Width`: integer pixels, default 512.
- `Color Table`: `CreateColorMapParameter`, default `Rainbow`.
- `Scale Mode`: `Automatic` or `Manual`.
- `Manual Minimum MUD`: active only in Manual mode.
- `Manual Maximum MUD`: active only in Manual mode.

#### Output PNG

- `Write PNG`: boolean, default true.
- `Output PNG File`: output-file path, active when `Write PNG` is true.

#### Output Image Geometry

- `Save as Image Geometry`: boolean, default true.
- `Output Image Geometry`: creation path, active when geometry output is true.
- The output contains a UInt8, three-component RGB array under `Cell Data`.

At least one output mode must be enabled.

### Preflight

Preflight validates:

- ODF geometry type, origin, uniform spacing, dimensions, and full-cube extent.
- Scalar Float32 or Float64 ODF array placement and tuple count.
- Equal tuple counts for Crystal Structures and Material Names.
- Phase index range and a supported EbsdLib crystal-structure value.
- Section count of at least 2.
- Sections per row from 1 through the section count.
- A section width of at least 128 pixels.
- A valid RGB color-table preset with at least two control colors.
- Finite manual limits with `0 <= minimum < maximum`.
- At least one enabled output.
- A `.png` destination when disk output is enabled.

Preflight reports:

- Selected phase, material, and Laue group.
- Input and rendered units.
- Full grid dimensions and spacing.
- Euler plotting limits.
- Exact section angles.
- Section layout and final composite dimensions.
- Automatic or manual scale mode.

When Image Geometry output is enabled, preflight creates the correctly sized
geometry and RGB array. It does not need the data-dependent automatic maximum to
determine image dimensions.

### Execution

Execution:

1. Checks cancellation.
2. Converts Float32 or Float64 input to an EbsdLib double buffer.
3. Resolves and validates the selected color-table control points.
4. Builds the EbsdLib configuration and calls `ODFSectionCompositor`.
5. Maps EbsdLib diagnostics to contextual SIMPLNX errors.
6. Copies RGBA to the optional three-component RGB Image Geometry.
7. Writes the optional PNG through `PngWriter` and `AtomicFile`.

Cancellation is checked between each stage. No destination file is replaced
until PNG encoding succeeds and the atomic file commits.

## Error behavior

- Negative or non-finite ODF values are errors; they are never silently
  clamped.
- An all-zero ODF is an execution error in automatic scale mode.
- Unknown or unsupported crystal-structure values are errors.
- Cropped, nonuniform, transposed, or incomplete Euler grids are rejected with
  actual dimensions, origin, spacing, and expected values in the message.
- Output file errors include the complete requested path and underlying writer
  diagnostic.
- An EbsdLib error includes its code and message plus the selected ODF path and
  phase context.

## Testing

### EbsdLib tests

- Hand-derived Count-Density-to-MUD conversion and SO(3) normalization.
- MUD passthrough.
- Exact section-angle generation for every supported Laue class.
- Expected rectangular Euler plotting limits for every Laue class.
- Exact linear interpolation between cell-centered phi2 planes.
- Exact periodic interpolation across 360 degrees to 0 degrees.
- Automatic and manual scale behavior.
- Endpoint color clipping and RGB control-point interpolation.
- Deterministic layout dimensions and applied metadata.
- Selected raw RGBA pixels covering panel data, background, axes, text regions,
  and color-bar endpoints.
- Invalid dimensions, spacing, units, symmetry, values, scale, section count,
  row count, and image size.

Tests use small analytical grids so the expected section values are independent
of the implementation.

### MTRSim tests

- Preflight actions and updated values.
- Full parameter-validation and error-path coverage.
- Float32 and Float64 input paths.
- Mathematically equivalent MUD and Count-Density inputs produce identical RGB
  output.
- Disk-only, geometry-only, and combined output modes.
- PNG dimensions and decoded pixels match the RGB Image Geometry.
- Atomic file behavior preserves an existing destination on failure.
- End-to-end `ComputeODFFilter` to `WriteODFSectionsFilter` on a small,
  deterministic HCP dataset.

### Validation dataset

Use `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1.ang` for repository-external
validation. Retain the generated EbsdLib/MTRSim montage, the MTEX 6.1.0 montage,
the pipeline, logs, and comparison notes outside the repository. The 324 MB ANG
file is not added to unit-test data.

The validation checks:

- Matching indexed-point population.
- Expected section angles and plotting limits.
- MUD color-scale values.
- Visual agreement of feature locations with MTEX.
- PNG and Image Geometry pixel agreement for the MTRSim output.

## Documentation

- Add the `WriteODFSectionsFilter` user-facing filter document.
- Add a canonical montage generated from the validation dataset.
- Document axis mapping, full-grid requirements, unit conversion, phase
  provenance, symmetry plotting regions, section interpolation, color scaling,
  and both output modes.
- Cross-link the new filter from `ComputeODFFilter` documentation.
- Add Doxygen documentation for the new public EbsdLib configuration, result,
  utilities, and compositor.

## Cross-repository delivery sequence

1. Implement and test the EbsdLib utilities and public compositor on the
   EbsdLib feature branch.
2. Implement and test the MTRSim filter against the local EbsdLib source using
   DREAM3D-NX preset `NX-Com-Qt69-Vtk96-Rel-EbsdLib`.
3. Run the external `1.ang` validation and review the generated montage against
   MTEX.
4. Once an EbsdLib release containing the public API is available, verify the
   normal `NX-Com-Qt69-Vtk96-Rel` dependency path.

No simplnx OrientationAnalysis source changes are required for v1.

# Compute ODF From Euler Angles

## Group (Subgroup)

MTRSim (Compute)

## Description

This **Filter** builds a discretized Orientation Distribution Function (ODF) on a regular Bunge-Euler-space grid from per-voxel EBSD data. It mirrors the algorithm in the MTRSim research code's MATLAB implementation (`matlab/calc_ODF.m`) and produces a single `Float64`, single-component **Cell** **Data Array** on a (new or existing) ODF **Image Geometry**.

There are two output modes:

- **Create New ODF Geometry** -- a brand-new **Image Geometry** is created from the user-supplied bin size, and the new ODF component is written into its cell **Attribute Matrix**.
- **Append to Existing ODF Geometry** -- a new component is added to the cell **Attribute Matrix** of an already-existing ODF **Image Geometry**. The bin size and grid dimensions are derived from that geometry's spacing; the user's `Bin Size [deg]` input is ignored. Use this to compute multiple ODFs (multi-region or multi-phase) onto the same Euler-space grid for direct comparison.

### Axis Mapping (Important)

The **Image Geometry** dimensions follow the same convention as the Read / Write MTRSim ODF (HDF5) **Filters**:

| Image axis | Bunge angle | Role            |
|------------|-------------|-----------------|
| **X**      | `phi2`      | fastest-varying |
| **Y**      | `PHI`       | middle          |
| **Z**      | `phi1`      | slowest-varying |

Spacing on every axis is in **degrees** and equal to the bin width. With the default 5 deg bin size, the grid is `(X, Y, Z) = (72, 36, 72)` for a total of `186,624` bins.

Standard ODF plots ("phi1-PHI plane at fixed phi2") correspond to the **YZ image-slicing plane** in DREAM3D-NX.

## Algorithm

For every voxel `i` where `phases[i] != 0` and (mask is unused OR `mask[i] == true`):

1. Read the Bunge Euler triple `(phi1, PHI, phi2)` (radians).
2. Look up the crystal structure with `crystal_structures[phases[i]]` (an EbsdLib `CrystalStructure::*` code, e.g. `Hexagonal_High = 0`, `Cubic_High = 1`).
3. Expand to all symmetrically-equivalent Bunge tuples using EbsdLib symmetry operators (12 for HCP, 24 for cubic).
4. Each symmetric tuple is converted to integer bin indices and deposited into the ODF grid:
   - **Smoothing off:** the entire deposit (1.0) goes into the center bin.
   - **Smoothing on:** the deposit is distributed across 27 bins using the tri-linear face/edge/corner weights from `calc_ODF.m`:
     - `0.332` to the center bin
     - `0.448` total spread evenly across the 6 face neighbors
     - `0.16`  total spread evenly across the 12 edge neighbors
     - `0.06`  total spread evenly across the 8 corner neighbors

After all voxels have been processed, the ODF is normalized by the **total count of symmetric-equivalent orientation deposits** (i.e. for each contributing voxel, the number of symmetric variants its phase produces; e.g. 12 per voxel for HCP, 24 per voxel for cubic). This matches the MATLAB `calc_ODF.m` convention (line 80: `N = size(phi1_vec, 1)` after symmetric expansion).

### Output Units (Important)

The **Filter** offers two output modes via the **Output Units** parameter. Both branches go through the same accumulator and symmetric-expansion math; the only difference is whether the per-bin Bunge-volume Jacobian is applied at the end.

**Count-Density** (default) -- the raw normalized histogram. Sums to `1.0` over all bins. Each deposit is divided by the total number of symmetric-equivalent orientation deposits across all contributing voxels. Bit-exact match to MATLAB `calc_ODF.m`. Use this whenever the goal is parity with the MATLAB reference, or when downstream consumers expect a Bunge-grid count distribution.

**MUD** (Multiples of Uniform Density) -- the SO(3) probability density divided by the uniform-orientation density. Use this for visualization that should look like MTEX's `plot(odf)`: sub-uniform regions read `< 1`, peak texture reads several MUD. After the count-density normalization, each `(i, j, k)` bin is rescaled by the per-row Jacobian factor so the result is comparable across the full Φ range:

```
MUD[i, j, k] = CountDensity[i, j, k] * 8 * pi^2 / (step_rad^3 * sin(PHI_center(j)))
```

where `step_rad = binSizeDeg * pi / 180` and `PHI_center(j) = (j + 0.5) * step_rad`. The factor amplifies low-Φ bins (which cover much less SO(3) volume) and shrinks equatorial bins, removing the gimbal-pole compression that makes a count-density plot look "empty" near `PHI = 0` and `PHI = pi` even when those orientations are physically present. For HCP α-Ti scans this typically converts a "blue at the top of the φ₁-Φ panel" Count-Density rendering into the MTEX-style "peaks visible at all Φ" view.

Without this Jacobian, a global linear colormap suppresses sub-peak features because Bunge bins near the Φ poles cover orders of magnitude less orientation space than equatorial bins -- equal counts represent very different physical densities.

### Smoothing Boundary Handling

When smoothing is enabled, all three Bunge axes use uniform modulo wrap for the +1/-1 smoothing-neighbor stencil -- matching `calc_ODF.m` lines 95-113. Bin **assignment** itself uses clamp-to-last-bin at the upper boundary (a sample exactly at `PHI = pi`, `phi1 = 2 pi`, or `phi2 = 2 pi` is placed in the LAST bin on that axis, not wrap-around to bin 0) -- matching `calc_ODF.m` lines 43-48. Validated bin-by-bin against the MATLAB reference output.

### Parallelization

Per-voxel work runs under SIMPLNX's `ParallelDataAlgorithm`. Each thread accumulates into its own private grid and the per-thread grids are merged after the parallel loop. There is no shared writable state during the parallel phase.

## Output

A single `Float64`, single-component **Cell** **Data Array** of size `nphi1 * nPHI * nphi2` (e.g., `72 * 36 * 72 = 186,624` at 5 deg spacing) on the target ODF **Image Geometry**, named per the **Component Name** parameter (default `Component 1`).

In **Create New** mode the **Filter** also creates the **Image Geometry** and the cell **Attribute Matrix** at the user-supplied paths. In **Append** mode the existing geometry's cell **Attribute Matrix** is reused, and the proposed component name must not collide with an existing array on it.

## Cropping the ODF to its Fundamental Zone

The output of `ComputeODFFilter` always spans the full Bunge-Euler cube
(`phi1 in [0 deg, 360 deg)`, `PHI in [0 deg, 180 deg]`, `phi2 in [0 deg, 360 deg)`;
dimensions 72 × 36 × 72 at the default 5 deg spacing). Most of this cube is
redundant under crystal symmetry — the same density information is replicated
symmetrically across multiple regions of Euler space. To produce a
fundamental-zone-only visualization comparable to MTEX's `plot(odf)` rendering,
follow `ComputeODFFilter` with the **Crop Image Geometry** filter using the
inclusive zero-based index ranges in the table below.

### Axis mapping reminder

| ImageGeom axis | Euler component | Full-cube range    | Full-cube cells at 5 deg |
|----------------|-----------------|--------------------|--------------------------|
| **X**          | `phi2`          | `[0 deg, 360 deg)` | 72                       |
| **Y**          | `PHI`           | `[0 deg, 180 deg]` | 36                       |
| **Z**          | `phi1`          | `[0 deg, 360 deg)` | 72                       |

`phi1` is never reduced by symmetry, so the Z axis is always cropped 0 - 71
(i.e., not cropped at all).

### Crop ranges by Laue class (5 deg spacing)

The EbsdLib code is the integer stored in the per-phase `CrystalStructures`
ensemble array (the same value `ComputeODFFilter` consumes). Crop bounds are
inclusive zero-based indices.

| EbsdLib Laue Class                 | Laue class (HM symbol) | X min / max (`phi2`) | Y min / max (`PHI`) | Z min / max (`phi1`) | Cropped dims (X x Y x Z) |
|------------------------------------|------------------------|----------------------|---------------------|----------------------|--------------------------|
| `0` — Hexagonal_High               | `6/mmm`                | 0 / 11               | 0 / 17              | 0 / 71               | 12 x 18 x 72             |
| `1` — Cubic_High                   | `m-3m`                 | 0 / 17               | 0 / 17              | 0 / 71               | 18 x 18 x 72             |
| `2` — Hexagonal_Low                | `6/m`                  | 0 / 11               | 0 / 35              | 0 / 71               | 12 x 36 x 72             |
| `3` — Cubic_Low                    | `m-3`                  | 0 / 35               | 0 / 17              | 0 / 71               | 36 x 18 x 72             |
| `4` — Triclinic                    | `-1`                   | 0 / 71               | 0 / 35              | 0 / 71               | 72 x 36 x 72 *(no crop)* |
| `5` — Monoclinic (`2 ‖ c` setting) | `2/m`                  | 0 / 35               | 0 / 35              | 0 / 71               | 36 x 36 x 72             |
| `5` — Monoclinic (`2 ‖ b` setting) | `2/m`                  | 0 / 71               | 0 / 17              | 0 / 71               | 72 x 18 x 72             |
| `6` — OrthoRhombic                 | `mmm`                  | 0 / 35               | 0 / 17              | 0 / 71               | 36 x 18 x 72             |
| `7` — Tetragonal_Low               | `4/m`                  | 0 / 17               | 0 / 35              | 0 / 71               | 18 x 36 x 72             |
| `8` — Tetragonal_High              | `4/mmm`                | 0 / 17               | 0 / 17              | 0 / 71               | 18 x 18 x 72             |
| `9` — Trigonal_Low                 | `-3`                   | 0 / 23               | 0 / 35              | 0 / 71               | 24 x 36 x 72             |
| `10` — Trigonal_High               | `-3m`                  | 0 / 23               | 0 / 17              | 0 / 71               | 24 x 18 x 72             |

### Derivation

Each Laue class' Euler-space FZ is determined by two symmetries:

- **`PHI` extent** is `[0 deg, 90 deg]` when the class has a 2-fold rotation
  perpendicular to the c-axis (i.e., the `mmm`, `m-3`, `m-3m`, `4/mmm`,
  `6/mmm`, `-3m` classes), and `[0 deg, 180 deg]` otherwise.
- **`phi2` extent** is `360 deg / N` where N is the order of the rotation
  about the c-axis: N = 6 (hex), 4 (tetragonal), 3 (trigonal), 2
  (orthorhombic and cubic-low along c, monoclinic in `2 ‖ c` setting),
  1 (triclinic and monoclinic in `2 ‖ b` setting).

### Monoclinic setting

Two crop options are listed for monoclinic because the FZ extents depend
on which crystallographic axis carries the 2-fold rotation:

- **`2 ‖ c` setting** — the 2-fold rotation is along the c-axis. The c-axis
  2-fold acts directly on `phi2`, halving its extent to `[0 deg, 180 deg)`.
  `PHI` keeps its full `[0 deg, 180 deg]` extent because no perpendicular
  2-fold is present. Used by some texture-analysis references that adopt
  c as the unique axis.
- **`2 ‖ b` setting** — the standard crystallographic setting in which b is
  the unique axis perpendicular to c. The b-axis 2-fold acts
  perpendicularly to the c-axis and halves `PHI` to `[0 deg, 90 deg]`;
  `phi2` keeps its full triclinic extent because there is no rotation
  symmetry along c.

Both rows describe the same point group; the choice depends on how the
monoclinic phase was set up upstream (e.g., in MTEX, in the EBSD vendor's
acquisition software, or in the `CreateEnsembleInfoFilter`). Pick the row
matching your dataset's setting. EbsdLib's `MonoclinicOps` uses one
specific setting internally; if you don't know which, check the operator
table in `EbsdLib/LaueOps/MonoclinicOps.cpp` or test by running
`ComputeODFFilter` on a known single-orientation input and observing
where the symmetric variants land.

### Cubic note

Cubic Laue classes (`m-3`, `m-3m`) are handled as rectangular bounding
boxes around the true tetrahedral cubic FZ. The bounding box is
approximately 3x larger than the strict FZ for `m-3m`, so the cropped
output still contains some rotational redundancy from the cubic 3-folds
along `<111>`. Eliminating that remaining redundancy would require a
non-rectangular crop, which is outside the scope of the SIMPLNX
**Crop Image Geometry** filter; the rectangular bounding box is the
standard convention used by MTEX, Bunge, and most ODF-visualization
tools.

### Multiphase ODF caveat

If the input scan contains multiple phases with different Laue classes
(e.g., an HCP `alpha` + BCC `beta` titanium alloy), no single rectangular
crop is optimal for both phases. In that case either compute one ODF per
phase (writing each into a separate ImageGeom) and crop each to its own
FZ, or leave the multiphase ODF uncropped.

### Preflight Preview

The UI shows the following before the user clicks Apply:

- Output Mode
- Target ODF Image Geometry
- Bin Size [deg] (and a note that it was derived from the existing geometry, in Append mode)
- Bins (phi1 x PHI x phi2)
- ImageGeom Dimensions (X, Y, Z)
- Total Bin Count
- Smoothing (enabled / disabled)
- Output Units (Count-Density / MUD)
- Mask (none / path)


% Auto generated parameter table will be inserted here

## Reference

`calc_ODF.m`, `symmetric_euler_angles.m` -- MATLAB sources distributed under `matlab/` in this plugin's repository (the original MTRSim research code). This **Filter** reproduces their grid layout, smoothing weights, and normalization convention.

## Example Pipelines

(Example pipelines will be added in a future release.)

## License & Copyright

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists

If you need help, need to file a bug report or want to request a new feature, please head over to the [DREAM3DNX-Issues](https://github.com/BlueQuartzSoftware/DREAM3DNX-Issues/discussions) GitHub site where the community of DREAM3D-NX users can help answer your questions.

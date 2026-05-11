# Write MTRSim ODF (HDF5)

## Group (Subgroup)

MTRSim (IO/Output)

## Description

This **Filter** writes a SIMPLNX **Image Geometry** plus a user-selected list of `Float64`, single-component **Cell** **Data Arrays** to a MATLAB-format MTRSim Orientation Distribution Function (ODF) HDF5 file. The output file is the same on-disk format produced by the MTRSim research code's MATLAB tooling and is consumed by the Read MTRSim ODF (HDF5) **Filter**, so a Read -> Write -> Read cycle is lossless.

### Axis Mapping (Important)

The on-disk file is keyed by Bunge-Euler axes (`phi1`, `PHI`, `phi2`), with `phi1` as the slowest-varying axis. The selected SIMPLNX **Image Geometry** is assumed to follow the same convention used by the Read **Filter**:

| Image axis | Bunge angle | Role on disk     |
| ---------- | ----------- | ---------------- |
| **X**      | `phi2`      | fastest-varying  |
| **Y**      | `PHI`       | middle           |
| **Z**      | `phi1`      | slowest-varying  |

Spacing on the **Image Geometry** is interpreted in **degrees**. The on-disk bin arrays (`phi1_bins`, `PHI_bins`, `phi2_bins`) are derived from the geometry's spacing on each axis. Spacing should be uniform across all three axes (this is not strictly enforced, but the file format assumes it).

### Output File Layout

For each selected component array, the **Filter** writes:

```
<prefix>/num_components                     (scalar, total component count)
<prefix>/component_0/ODFval                 (flattened Float64 ODF values)
<prefix>/component_0/phi1_bins              (Float64 bin centers, degrees)
<prefix>/component_0/PHI_bins               (Float64 bin centers, degrees)
<prefix>/component_0/phi2_bins              (Float64 bin centers, degrees)
<prefix>/component_1/...
...
```

Default prefix is `/ODF_best` to match the MATLAB convention.

### Selection Order Determines Component Numbering (Important)

The order in which arrays appear in the **ODF Component Arrays** selection is the **source of truth** for the file's `component_N` numbering. The first array in the selection becomes `component_0` on disk, the second becomes `component_1`, and so on. The original array names in the **DataStructure** are not preserved in the output file.

For example, if you import a 3-component ODF, then later select `component_2` and `component_5` from the imported geometry in that order, the exported file's `component_0` will hold the input's `component_2` data and `component_1` will hold the input's `component_5` data. Re-order the selection in the UI if you want a different mapping.

### Preflight Preview

Preflight populates the following so the UI can preview what will be written:

- HDF5 Path Prefix
- Components to Write
- Tuples per Component
- Grid (phi1 x PHI x phi2)
- Spacing (phi1, PHI, phi2) [deg]
- Estimated File Size [bytes]

## Notes

- An existing file at the output path is overwritten without warning.
- The selected geometry must be an **Image Geometry**; other geometry types are rejected.
- Each selected array must live on the cell **Attribute Matrix** of the selected geometry, must be `Float64`, must be single-component, and must have a tuple count equal to the geometry's cell count.
- All three image-axis dimensions must be at least 2 (degenerate single-slice grids are rejected).

## Errors

| Code     | Meaning                                                                                  |
| -------- | ---------------------------------------------------------------------------------------- |
| `-12100` | No component arrays selected.                                                            |
| `-12101` | Selected geometry is not an **Image Geometry**.                                          |
| `-12102` | A selected array does not belong to the selected **Image Geometry**.                     |
| `-12103` | A selected array is not a `Float64` **Data Array**.                                      |
| `-12104` | A selected array's tuple count does not match the geometry's cell count.                 |
| `-12105` | Degenerate geometry (fewer than 2 bins on at least one axis).                            |
| `-12110` | The HDF5 write call failed; error message includes the underlying cause.                 |

% Auto generated parameter table will be inserted here

## Example Pipelines

(Example pipelines will be added in a future release.)

## License & Copyright

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists

If you need help, need to file a bug report or want to request a new feature, please head over to the [DREAM3DNX-Issues](https://github.com/BlueQuartzSoftware/DREAM3DNX-Issues/discussions) GitHub site where the community of DREAM3D-NX users can help answer your questions.

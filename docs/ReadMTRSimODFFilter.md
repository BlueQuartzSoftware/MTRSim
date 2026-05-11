# Read MTRSim ODF (HDF5)

## Group (Subgroup)

MTRSim (IO/Input)

## Description

This **Filter** reads a MATLAB-format MTRSim Orientation Distribution Function (ODF) HDF5 file into the **DataStructure**. It creates a single **Image Geometry** representing the discretized Bunge-Euler space and one `Float64`, single-component **Cell** **Data Array** per ODF component stored in the file.

The on-disk layout matches the convention used by the MTRSim research code's MATLAB tooling (see `matlab/calc_ODF.m` and the companion writer routines). A file written by the Write MTRSim ODF (HDF5) **Filter** can be round-tripped losslessly with this **Filter**.

### Axis Mapping (Important)

The file stores the ODF on a regular grid in Bunge-Euler space (`phi1`, `PHI`, `phi2`), with `phi1` as the slowest-varying axis and `phi2` as the fastest. To preserve the row-major in-memory ordering when materializing the data into a SIMPLNX **Image Geometry** (whose tuple order is Z, Y, X), this **Filter** maps the Euler axes to image axes as follows:

| Bunge angle | Role on disk          | Image axis |
| ----------- | --------------------- | ---------- |
| `phi2`      | fastest-varying       | **X**      |
| `PHI`       | middle                | **Y**      |
| `phi1`      | slowest-varying       | **Z**      |

The **Image Geometry**'s spacing is set in **degrees** on every axis and is equal to the bin width.

**UX consequence:** the standard ODF plot ("phi1-PHI plane at fixed phi2") corresponds to DREAM3D-NX's **YZ image-slicing plane** at a chosen X slice. If you are used to seeing phi1-PHI sections in MATLAB, slice along X in the visualization, not along Z.

### Outputs

- An **Image Geometry** at the user-selected path with dimensions `(X, Y, Z) = (phi2 bins, PHI bins, phi1 bins)` and uniform spacing in degrees.
- A **Cell** **Attribute Matrix** at the user-selected name under the **Image Geometry**.
- One `Float64`, single-component **Data Array** per component in the file, named `component_0`, `component_1`, ..., `component_{N-1}`.

### Preflight Preview

The **Filter** opens the file during preflight to read its metadata, so the UI shows the following before the user clicks Apply:

- HDF5 Path Prefix
- Components Found
- Bins (phi1 x PHI x phi2)
- ImageGeom Dimensions (X, Y, Z)
- Spacing (X, Y, Z) [deg]

Preflight also validates that:

- The file exists and is a readable HDF5 file.
- `<prefix>/num_components` is at least 1.
- All `N` components share byte-exact identical bin arrays (`phi1_bins`, `PHI_bins`, `phi2_bins`).
- The bin arrays are uniformly spaced and monotonically increasing.
- The `ODFval` array sizes match the dimensions derived from the bin arrays.

## Notes

- The default HDF5 path prefix is `/ODF_best`, matching the variable name MATLAB workflows typically use. If your file was written under a different MATLAB variable name, set the prefix to `/<that_name>`.
- All ODF arrays are read as `Float64` regardless of how they are stored on disk.

## Errors

| Code     | Meaning                                                                                                                                |
| -------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `-12001` | (Preflight) File read failed. The error message includes the underlying cause (file missing, wrong format, mismatched bin arrays across components, non-uniform bins, dimension mismatch, etc.). |
| `-12010` | (Execute) File read failed after preflight succeeded (e.g., the file became unreadable between preflight and execute). |
| `-12012` | (Execute) Per-component array size in the file does not match the array created during preflight. |

% Auto generated parameter table will be inserted here

## Example Pipelines

(Example pipelines will be added in a future release.)

## License & Copyright

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists

If you need help, need to file a bug report or want to request a new feature, please head over to the [DREAM3DNX-Issues](https://github.com/BlueQuartzSoftware/DREAM3DNX-Issues/discussions) GitHub site where the community of DREAM3D-NX users can help answer your questions.

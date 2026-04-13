# MTRsim

Copyright (c) 2025 US Govt

MTRsim generates synthetic volumes of microstructures containing **Microtexture Regions
(MTRs)**. Given a set of component Orientation Distribution Functions (ODFs) with
corresponding volume fractions and spatial correlation parameters, and the dimensions of
an output voxel grid, the tool produces a simulated random draw of a volume containing
spatially coherent regions of similar crystallographic orientation.

Distribution A. Approved For Public Release: Distribution Unlimited. AFRL-2025-1631

---

## Disclaimer

The software is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE WARRANTY OF
DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

---

## Feedback

Daniel.Sparkman.1 at us dot af dot mil

---

## Building

### Prerequisites

- CMake 3.21+
- Ninja
- vcpkg (path configured in `CMakeUserPresets.json`)
- EbsdLib built and installed at the prefix referenced by `CMakeUserPresets.json`

### Configure and build

```bash
# Configure (reads CMakeUserPresets.json for vcpkg and EbsdLib paths)
cd /path/to/MTR_sim
cmake --preset mtrsim-Rel

# Build the library and executable
cd /path/to/Build/mtrsim-Rel
cmake --build . --target all

# Run unit tests
ctest --test-dir .
```

The executable is written to `<build-dir>/EbsdLib/Bin/mtrsim`.

---

## Data Setup (one-time)

The simulation reads ODF data from `data/simulation_ODF.h5`. This file is produced by
converting the original MATLAB v5 `.mat` file using the provided Python script:

```bash
# Requires: Python with h5py and scipy (conda environment 'dream3d' works)
conda activate dream3d
python tools/convert_simulation_ODF.py
```

The remaining data files (`blank_ODF.mat`, `uniformODF.mat`, `PF_coords_for_ODF.mat`,
`PF_indexed.mat`) were already in MATLAB v7.3 (HDF5) format and do not require
conversion.

---

## Running the Executable

### Synopsis

```
mtrsim [OPTIONS]

OPTIONS:
  -h, --help              Print help and exit
  -c, --config <path>     Path to a JSON configuration file
  -o, --output <dir>      Output directory (default: current directory)
      --seed <uint>        Fixed RNG seed for reproducible runs (0 = random_device)
```

### Quickstart — all defaults

```bash
./mtrsim
```

Runs with the default parameters from `simulate_MTRs.m` (see table below), reading the
ODF from `data/simulation_ODF.h5` relative to the working directory, and writing outputs
to the current directory.

### Reproducible run

```bash
./mtrsim --seed 42
```

### Custom config and output directory

```bash
./mtrsim -c params.json -o results/
```

`--seed` on the command line takes precedence over a `seed` key in the JSON file.

---

## Configuration

All simulation parameters can be supplied via a JSON file passed to `-c/--config`. Every
key is optional; omitted keys use the defaults below.

### JSON keys and defaults

| Key | Type | Default | Description |
|---|---|---|---|
| `xLen` | float | `38.1` | Volume length in x [mm] (`1.5 * 25.4`) |
| `yLen` | float | `12.7` | Volume length in y [mm] (`0.5 * 25.4`) |
| `zLen` | float | `0.0` | Volume length in z [mm] (0 = 2-D simulation) |
| `dx` | float | `0.02` | Voxel spacing in x [mm] |
| `dy` | float | `0.02` | Voxel spacing in y [mm] |
| `dz` | float | `0.02` | Voxel spacing in z [mm] |
| `volumeFractions` | float[] | `[0.30, 0.35, 0.35]` | Volume fraction of each MTR component ODF (must sum to 1) |
| `thetaList` | float[][] | `[[0.10,0.45,0.10],[0.08,0.37,0.08]]` | Correlation length parameters — rows are latent Gaussian fields, columns are x/y/z directions. Increasing a value increases MTR size in that direction. |
| `nuggetVariance` | float[] | `[0.67, 0.71, 0.72]` | Nugget variance per component (currently unused in simulation) |
| `odfInputPath` | string | `"data/simulation_ODF.h5"` | Path to the HDF5 ODF input file |
| `seed` | uint | `0` | RNG seed (0 = non-deterministic via `std::random_device`) |

### Example — full config matching MATLAB defaults

```json
{
  "xLen": 38.1,
  "yLen": 12.7,
  "zLen": 0.0,
  "dx": 0.02,
  "dy": 0.02,
  "dz": 0.02,
  "volumeFractions": [0.30, 0.35, 0.35],
  "thetaList": [
    [0.10, 0.45, 0.10],
    [0.08, 0.37, 0.08]
  ],
  "nuggetVariance": [0.67, 0.71, 0.72],
  "odfInputPath": "data/simulation_ODF.h5",
  "seed": 0
}
```

### Grid dimensions

The voxel grid is derived from the parameters as:

```
nx = round(xLen / dx)
ny = round(yLen / dy)
nz = max(round(zLen / dz), 1)   # always at least 1 z-layer
N  = nx * ny * nz
```

With the defaults this gives `nx=1905 × ny=635 × nz=1 = 1,209,675` voxels.

---

## Outputs

Both output files are written to the directory specified by `-o` (default: `.`).

| File | Description |
|---|---|
| `sim_IPF_map.png` | Inverse Pole Figure (IPF) colour map image. HCP crystal symmetry, Z-axis reference direction. |
| `sim_results.csv` | Per-voxel data table (see schema below). |

### `sim_results.csv` schema

| Column | Units | Description |
|---|---|---|
| `x` | mm | Voxel x-coordinate |
| `y` | mm | Voxel y-coordinate |
| `z` | mm | Voxel z-coordinate |
| `phi1` | rad | Bunge Euler angle φ₁ |
| `PHI` | rad | Bunge Euler angle Φ |
| `phi2` | rad | Bunge Euler angle φ₂ |
| `mtr_index` | — | 1-based MTR component assignment |

Voxel ordering is **z outer → x middle → y inner**, matching `simulate_MTRs.m`:

```
row k  →  x = j·dx,  y = i·dy,  z = zix·dz
         j ∈ [1, nx],  i ∈ [1, ny],  zix ∈ [1, nz]
```

---

## Method Overview

The simulation proceeds in four stages:

1. **Plurigaussian Random Field (PGRF)** — two correlated Gaussian random fields are
   generated over the voxel grid using a Kronecker-structured Cholesky decomposition
   (`gp_generator.m` → `GPGenerator`). Thresholds derived from the volume fractions
   assign each voxel to one of the MTR component classes (`select_AR.m` / `eval_AR.m`
   → `AssignmentRule`).

2. **ODF loading** — component ODFs are read from `simulation_ODF.h5`. Each ODF is a
   discrete histogram over the 72×36×72 Bunge Euler-angle grid (5° spacing,
   186,624 bins). ODFs are normalised to sum to 1.

3. **Orientation sampling** — for each component, N orientations are drawn by
   inverse-CDF sampling of the ODF histogram (`sample_N_orientations_from_ODF.m`
   → `ODFSampler`). Each voxel receives the orientation sampled for its assigned
   component.

4. **Output** — an IPF colour map PNG is rendered (`view_IPF_map.m` → `IPFMapper`) and
   the full per-voxel orientation table is written to CSV.

---

## Source File Mapping

| MATLAB file | C++ module | Role |
|---|---|---|
| `simulate_MTRs.m` | `app/main.cpp` | Top-level driver |
| `PGRF_simulation.m` | `PGRFSimulation` | Plurigaussian RF simulation |
| `gp_generator.m` | `GPGenerator` | Gaussian process field via Kronecker Cholesky |
| `select_AR.m` / `eval_AR.m` | `AssignmentRule` | Threshold selection and voxel classification |
| `qsimvn.m` | `QSimVN` | Quasi-Monte Carlo multivariate normal CDF |
| `calc_ODF.m` | `ODFCalculator` | ODF histogram from Euler angles |
| `sample_orientation_from_ODF.m` / `sample_N_orientations_from_ODF.m` | `ODFSampler` | Draw orientations from ODF |
| `symmetric_euler_angles.m` | `CrystalSymmetry` | Symmetrically equivalent orientations (wraps EbsdLib) |
| `unit_triangle_IPF_coords.m` / `IPF_colors.m` / `view_IPF_map.m` | `IPFMapper` | IPF projection and PNG output |
| `convert_ODF_to_PF.m` | `PoleFigure` | ODF → pole figure intensities |
| `load_MTR_data.m` | `MTRDataLoader` | Load experimental EBSD data from CSV |

---

## Acknowledgements

The `qsimvn.m` provided is based on the function of the same name by **Alan Genz**. It
has been modified slightly for use with this software. The full acknowledgement, original
redistribution conditions, copyright notice, and disclaimer are contained in
`matlab/qsimvn.m`. In accordance with the redistribution conditions, no claim of
endorsement by Alan Genz is made for MTRsim.

---

## Citing

The software has not been published yet. The appropriate citation will be:

```bibtex
@article{MTRsim,
  title   = {Spatial Statistical Simulation of Microtexture},
  author  = {Daniel Sparkman, Laura Homa, Adam Pilchak, and Harry Millwater},
  journal = {###},
  volume  = {#},
  number  = {#},
  month   = {###.},
  year    = {2024}
}
```

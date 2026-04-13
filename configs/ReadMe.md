# MTRSim Configuration Files

This directory contains JSON configuration files for the MTRSim microtexture region simulator. Each file can be passed to both the C++ `mtrsim` executable and the MATLAB wrapper `run_simulation_from_json.m` to produce comparable outputs.

## Configuration Schema

All fields are optional. Omitted fields use the defaults shown below (matching the original MATLAB `simulate_MTRs.m` parameters).

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `xLen` | float | 38.1 | mm | Physical volume extent in x direction |
| `yLen` | float | 12.7 | mm | Physical volume extent in y direction |
| `zLen` | float | 0.0 | mm | Physical volume extent in z direction (0 = 2D simulation) |
| `dx` | float | 0.02 | mm | Voxel spacing in x |
| `dy` | float | 0.02 | mm | Voxel spacing in y |
| `dz` | float | 0.02 | mm | Voxel spacing in z |
| `volumeFractions` | float[] | [0.30, 0.35, 0.35] | — | Target volume fraction per MTR component (must sum to 1) |
| `thetaList` | float[][] | [[0.10,0.45,0.10],[0.08,0.37,0.08]] | mm | Correlation lengths for each latent Gaussian field. Rows = Gaussians (count = `len(volumeFractions) - 1`), columns = [theta_x, theta_y, theta_z] |
| `nuggetVariance` | float[] | [0.67, 0.71, 0.72] | — | Nugget variance per component (currently unused) |
| `odfInputPath` | string | "data/simulation_ODF.h5" | — | Path to HDF5 ODF data file |
| `seed` | uint64 | 0 | — | Random seed (0 = non-deterministic) |

**Key relationships:**
- Image dimensions: `nx = round(xLen/dx)`, `ny = round(yLen/dy)`, `nz = max(round(zLen/dz), 1)`
- Total voxels: `N = nx * ny * nz`
- Number of latent Gaussians: `len(volumeFractions) - 1`
- `thetaList` must have exactly `len(volumeFractions) - 1` rows

## Existing Configurations

### Production / Demonstration Configs

These use the full default domain (38.1 mm x 12.7 mm) and are suitable for generating publication-quality outputs.

| File | Image Size | Description |
|------|-----------|-------------|
| `default.json` | 1905 x 635 px | Baseline configuration matching the original MATLAB `simulate_MTRs.m` defaults. Standard needle-like MTR morphology with 3 components at 30%/35%/35% volume fractions. |
| `image_low_res.json` | 953 x 318 px | Same physical volume at half resolution (dx=dy=0.04 mm). ~8x faster. Useful for quick iteration. |
| `image_high_res.json` | 3810 x 1270 px | Same physical volume at double resolution (dx=dy=0.01 mm). 4x the pixels, significantly slower due to O(nx^2) Cholesky cost. |
| `image_large_volume.json` | 3810 x 1270 px | Double physical size (76.2 x 25.4 mm) at default resolution. Shows more MTR regions across a larger scan area. |
| `mtrs_large.json` | 1905 x 635 px | Theta values doubled — MTR regions are physically ~2x larger (theta_y ~ 0.90 mm). |
| `mtrs_small.json` | 1905 x 635 px | Theta values at ~40% of default — fine-grained texture with many small MTRs (theta_y ~ 0.18 mm). |
| `mtrs_isotropic.json` | 1905 x 635 px | Equal theta in all directions (0.30/0.25 mm). Produces blocky, equiaxed MTR regions instead of elongated needles. |
| `mtrs_very_elongated.json` | 1905 x 635 px | Extreme y-direction correlation (theta_y = 1.50 mm, theta_x = 0.05 mm). ~30:1 aspect ratio needle-like MTRs. |

### Test Configurations

These use small grids (typically 100x50 pixels or smaller) for fast automated testing. Designed to run in seconds rather than minutes.

| File | Image Size | What It Tests |
|------|-----------|---------------|
| `test_tiny_default.json` | 100 x 50 px | **Smoke test.** Same parameters as `default.json` on a tiny grid. Verifies the full pipeline runs end-to-end. |
| `test_equal_fractions.json` | 100 x 50 px | **Symmetric assignment rule.** Equal volume fractions (1/3 each). All components should appear in roughly equal proportions. |
| `test_dominant_component.json` | 100 x 50 px | **Asymmetric fractions.** 70%/15%/15% split. Component 1 forms large connected regions; components 2 and 3 appear as small islands. |
| `test_extreme_dominant.json` | 100 x 50 px | **Near-degenerate fractions.** 90%/5%/5% split. Stress-tests assignment rule thresholds pushed to extremes. |
| `test_3d_thin_slab.json` | 100 x 50 x 5 | **3D simulation.** Non-zero zLen (0.10 mm, 5 z-layers). Tests z-direction Cholesky decomposition and z-major voxel ordering. |
| `test_square_domain.json` | 100 x 100 px | **Square domain.** Equal xLen and yLen (2.0 mm each). Anisotropic theta should still produce elongated MTRs despite the symmetric domain. |
| `test_seed_99.json` | 100 x 50 px | **Seed sensitivity.** Same as `test_tiny_default.json` but seed=99. Results should differ spatially but have the same statistical properties. |
| `test_symmetric_theta.json` | 100 x 50 px | **Identical Gaussians.** Both latent fields use the same theta values. Removes inter-Gaussian morphological variation. |
| `test_high_aspect_theta.json` | 100 x 50 px | **Orthogonal elongation.** Gaussian 1 elongated in y, Gaussian 2 elongated in x. Creates complex cross-hatched MTR morphology. |
| `test_coarse_voxels.json` | 20 x 10 px | **Under-resolved regime.** Voxel size (0.10 mm) comparable to theta — only 200 voxels. Tests behavior when spatial resolution is poor. |
| `test_large_theta.json` | 100 x 50 px | **Domain-spanning correlation.** theta_y (5.0 mm) exceeds domain yLen (1.0 mm). Gaussian field nearly constant in y → stripe-like banding. |
| `test_small_theta.json` | 100 x 50 px | **White noise limit.** theta (0.001 mm) much smaller than voxel spacing. Covariance matrix approaches identity; PGRF degenerates to i.i.d. draws. Volume fractions should be correct but no spatial structure. |

## Running Simulations

### C++ (command line)

```bash
# Build first
cd /Users/mjackson/Workspace5/Build/mtrsim-Rel && cmake --build . --target all

# Run a single config
cd /path/to/MTRSim
./scripts/run_cpp.sh configs/test_tiny_default.json output/test_tiny_default

# Output files:
#   output/test_tiny_default/cpp_sim_results.csv
#   output/test_tiny_default/cpp_sim_IPF_map.png
```

### MATLAB (command line)

```bash
cd /path/to/MTRSim
./scripts/run_matlab.sh configs/test_tiny_default.json output/test_tiny_default

# Output files:
#   output/test_tiny_default/matlab_sim_results.csv
#   output/test_tiny_default/matlab_sim_IPF_map.jpg
#   output/test_tiny_default/matlab_sim_PGRF.csv
```

### Full Comparison Suite

```bash
# Run all test configs through both MATLAB and C++, then compare:
./scripts/run_comparison_suite.sh

# Run only MATLAB:
./scripts/run_comparison_suite.sh --matlab-only

# Run only C++:
./scripts/run_comparison_suite.sh --cpp-only

# Re-run comparison on existing outputs:
./scripts/run_comparison_suite.sh --compare-only

# Run a subset of configs:
./scripts/run_comparison_suite.sh --configs "configs/test_tiny*.json"
```

### Direct MATLAB Invocation

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -nodisplay -nosplash -batch \
    "addpath('matlab'); addpath('data'); run_simulation_from_json('configs/default.json', 'output/matlab')"
```

## Output Format

Both MATLAB and C++ produce CSV files with the same column layout:

```
x,y,z,phi1,PHI,phi2,mtr_index
0.020000,0.020000,0.020000,1.234567,0.567890,2.345678,1
0.020000,0.040000,0.020000,0.987654,0.123456,1.876543,2
...
```

| Column | Description |
|--------|-------------|
| `x, y, z` | Voxel spatial coordinates in mm |
| `phi1, PHI, phi2` | Bunge Euler angles in radians |
| `mtr_index` | Component assignment (1-based) |

Output files use prefixes to distinguish source:
- `matlab_sim_results.csv`, `matlab_sim_IPF_map.jpg` — MATLAB output
- `cpp_sim_results.csv`, `cpp_sim_IPF_map.png` — C++ output

## Comparison Notes

The MATLAB and C++ implementations use **different random number generators** (MATLAB's `mt19937ar` vs C++ `std::mt19937_64`). Even with the same seed, outputs will not be identical. The comparison script (`scripts/compare_results.py`) performs **statistical** validation:

1. **Grid geometry** — Same number of voxels, same coordinate ranges
2. **Component coverage** — All expected components present
3. **Volume fractions** — Component proportions match targets within 15% tolerance (correlated spatial fields produce more variance than i.i.d. samples, especially on small grids)
4. **Euler angle ranges** — All angles within valid bounds
5. **Euler angle distributions** — Per-component mean and standard deviation are in the same ballpark (informational, not a pass/fail criterion)

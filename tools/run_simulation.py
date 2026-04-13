#!/opt/local/anaconda3/envs/dream3d/bin/python3
"""
run_simulation.py — Python driver for MTRSim.

Replicates the full workflow of app/main.cpp using the mtrsim Python bindings:
  1. Configure simulation parameters
  2. Run the PGRF simulation (voxel-level MTR assignment)
  3. Load ODF components from HDF5
  4. Sample Euler-angle orientations per component
  5. Assign per-voxel orientations
  6. Write IPF-colour map PNG and results CSV

Usage:
    python tools/run_simulation.py [--odf PATH] [--output DIR] [--seed N] [--preset NAME]

Prerequisites:
    conda activate dream3d   (h5py and numpy are already present)
    export PYTHONPATH=/Users/mjackson/Workspace7/Build/mtrsim-Rel/src/libmtrsim

Example (from the MTRSim source directory):
    PYTHONPATH=/Users/mjackson/Workspace7/Build/mtrsim-Rel/src/libmtrsim \\
        /opt/local/anaconda3/envs/dream3d/bin/python3 tools/run_simulation.py \\
        --odf data/simulation_ODF.h5 --seed 42

    # Or, after activating the conda env and setting PYTHONPATH:
    conda activate dream3d
    export PYTHONPATH=/Users/mjackson/Workspace7/Build/mtrsim-Rel/src/libmtrsim
    python tools/run_simulation.py --preset default --seed 42
"""

import argparse
import csv
import math
import os
import sys

import h5py
import numpy as np

import mtrsim


# ---------------------------------------------------------------------------
# ODF helpers
# ---------------------------------------------------------------------------

def load_odf_components(hdf_path: str) -> list[mtrsim.ODFComponent]:
    """
    Load ODF components from an HDF5 file produced by convert_simulation_ODF.py.

    Expected layout:
        /ODF_best/num_components          — int64 scalar
        /ODF_best/component_N/ODFval      — (186624,) float64
        /ODF_best/component_N/phi1_bins   — (72,)     float64  [rad]
        /ODF_best/component_N/PHI_bins    — (36,)     float64  [rad]
        /ODF_best/component_N/phi2_bins   — (72,)     float64  [rad]

    Note: odfVal arrays are normalised so they sum to 1 before being passed to
    ODFSampler (matches the pre-processing in simulate_MTRs.m).
    """
    components = []
    with h5py.File(hdf_path, "r") as f:
        n = int(f["/ODF_best/num_components"][()])
        print(f"  ODF file: {n} component(s)")
        for j in range(n):
            pfx = f"/ODF_best/component_{j}"
            odf_val  = f[pfx + "/ODFval"][:]
            phi1     = f[pfx + "/phi1_bins"][:]
            phi      = f[pfx + "/PHI_bins"][:]
            phi2     = f[pfx + "/phi2_bins"][:]

            # Normalise so weights sum to 1
            total = odf_val.sum()
            if total > 0.0:
                odf_val = odf_val / total

            comp = mtrsim.ODFComponent()
            comp.odf_val   = odf_val
            comp.phi1_bins = phi1
            comp.phi_bins  = phi
            comp.phi2_bins = phi2
            components.append(comp)

    return components


def build_uniform_odf() -> mtrsim.ODFComponent:
    """
    Build the flat 186 624-element reference ODF used by ODFSampler.

    Bin centres follow the same formula as ODFCalculator::compute():
        phi1_bins[ix] = (i1   + 0.5) * 2π/72
        phi_bins[ix]  = (iPHI + 0.5) *  π/36
        phi2_bins[ix] = (i2   + 0.5) * 2π/72
    where ix = i1*(36*72) + iPHI*72 + i2.
    """
    n1, nPHI, n2 = 72, 36, 72
    n_total = n1 * nPHI * n2  # 186 624

    ix = np.arange(n_total)
    i1   =  ix // (nPHI * n2)
    iPHI = (ix  % (nPHI * n2)) // n2
    i2   =  ix  % n2

    uni = mtrsim.ODFComponent()
    uni.odf_val   = np.full(n_total, 1.0 / n_total)
    uni.phi1_bins = (i1   + 0.5) * (2.0 * math.pi / n1)
    uni.phi_bins  = (iPHI + 0.5) * (      math.pi / nPHI)
    uni.phi2_bins = (i2   + 0.5) * (2.0 * math.pi / n2)
    return uni


# ---------------------------------------------------------------------------
# Spatial coordinate grid
# ---------------------------------------------------------------------------

def build_spatial_coords(params: mtrsim.SimulationParams) -> tuple[int, int, int, np.ndarray]:
    """
    Build the [N x 3] spatial coordinate array (x, y, z) in mm.

    Ordering matches simulate_MTRs.m: z outer, x middle, y inner.
        s(k) = [j*dx, i*dy, zix*dz]
        j ∈ [1, nx],  i ∈ [1, ny],  zix ∈ [1, nz]
    """
    nx  = round(params.x_len / params.dx)
    ny  = round(params.y_len / params.dy)
    nz  = max(round(params.z_len / params.dz), 1)
    N   = nx * ny * nz

    coords = np.empty((N, 3), dtype=np.float64)
    k = 0
    for zix in range(1, nz + 1):
        for j in range(1, nx + 1):
            for i in range(1, ny + 1):
                coords[k, 0] = j   * params.dx
                coords[k, 1] = i   * params.dy
                coords[k, 2] = zix * params.dz
                k += 1

    return nx, ny, nz, coords


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_csv(path: str, spatial_coords: np.ndarray,
              phi1: np.ndarray, phi: np.ndarray, phi2: np.ndarray,
              mtr_index: np.ndarray) -> None:
    """Write per-voxel results to CSV (matches the format of app/main.cpp)."""
    N = len(phi1)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "y", "z", "phi1", "PHI", "phi2", "mtr_index"])
        for i in range(N):
            writer.writerow([
                f"{spatial_coords[i, 0]:.6f}",
                f"{spatial_coords[i, 1]:.6f}",
                f"{spatial_coords[i, 2]:.6f}",
                f"{phi1[i]:.6f}",
                f"{phi[i]:.6f}",
                f"{phi2[i]:.6f}",
                mtr_index[i],
            ])


# ---------------------------------------------------------------------------
# Simulation preset configurations
# ---------------------------------------------------------------------------

PRESETS = {
    "default": dict(
        x_len=38.1, y_len=12.7, z_len=0.0,
        dx=0.02,    dy=0.02,    dz=0.02,
        volume_fractions=[0.30, 0.35, 0.35],
        theta_list=[[0.10, 0.45, 0.10], [0.08, 0.37, 0.08]],
        nugget_variance=[0.67, 0.71, 0.72],
    ),
    "small": dict(
        x_len=5.0, y_len=5.0, z_len=0.0,
        dx=0.05,   dy=0.05,  dz=0.05,
        volume_fractions=[0.30, 0.35, 0.35],
        theta_list=[[0.10, 0.45, 0.10], [0.08, 0.37, 0.08]],
        nugget_variance=[0.67, 0.71, 0.72],
    ),
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="MTRSim Python driver")
    parser.add_argument("--odf",    default="data/simulation_ODF.h5",
                        help="Path to simulation_ODF.h5 (default: data/simulation_ODF.h5)")
    parser.add_argument("--output", default=".",
                        help="Output directory for PNG and CSV (default: .)")
    parser.add_argument("--seed",   type=int, default=42,
                        help="Random seed (0 = non-deterministic, default: 42)")
    parser.add_argument("--preset", choices=list(PRESETS), default="small",
                        help="Simulation size preset (default: small)")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # ── 1. Build simulation parameters ──────────────────────────────────────
    print("=== MTRSim Python driver ===")
    print(f"Preset : {args.preset}")
    print(f"ODF    : {args.odf}")
    print(f"Output : {args.output}")
    print(f"Seed   : {args.seed}")
    print()

    params = mtrsim.SimulationParams()
    preset = PRESETS[args.preset]
    params.x_len            = preset["x_len"]
    params.y_len            = preset["y_len"]
    params.z_len            = preset["z_len"]
    params.dx               = preset["dx"]
    params.dy               = preset["dy"]
    params.dz               = preset["dz"]
    params.volume_fractions = preset["volume_fractions"]
    params.theta_list       = preset["theta_list"]
    params.nugget_variance  = preset["nugget_variance"]
    params.odf_input_path   = args.odf
    params.output_dir       = args.output
    params.seed             = args.seed

    nx, ny, nz, spatial_coords = build_spatial_coords(params)
    N = nx * ny * nz
    print(f"Grid   : nx={nx}  ny={ny}  nz={nz}  N={N:,}")
    print()

    # ── 2. Run PGRF simulation ────────────────────────────────────────────────
    # Rng is a Python-owned wrapper around std::mt19937_64.  Pass it to every
    # stochastic class; py::keep_alive ensures the engine outlives all users.
    rng = mtrsim.Rng(args.seed)

    print("Running PGRF simulation...")
    pgrf   = mtrsim.PGRFSimulation(rng)
    result = pgrf.run(params)
    print(f"  mtr_index shape : {result.mtr_index.shape}")
    print(f"  latent_fields   : {result.latent_fields.shape}")

    # mtr_index is 1-based; count how many voxels ended up in each component
    unique, counts = np.unique(result.mtr_index, return_counts=True)
    for comp_id, cnt in zip(unique, counts):
        print(f"    component {comp_id}: {cnt:6d} voxels  ({100.0 * cnt / N:.1f} %)")
    print()

    # ── 3. Load ODF components from HDF5 ─────────────────────────────────────
    print(f"Loading ODF from: {args.odf}")
    odf_components = load_odf_components(args.odf)
    uniform_odf    = build_uniform_odf()
    n_components   = len(odf_components)
    print()

    # ── 4. Sample N orientations per component ────────────────────────────────
    # ODFSampler.sample_n() returns a [N x 3] numpy array: (phi1, PHI, phi2) [rad].
    sampler = mtrsim.ODFSampler(rng)
    orient_samples = []  # one [N x 3] array per component
    for j, comp in enumerate(odf_components):
        print(f"Sampling component {j + 1}/{n_components}  (N={N:,})...")
        samples = sampler.sample_n(N, comp, uniform_odf)
        orient_samples.append(samples)
    print()

    # ── 5. Assign per-voxel orientations ──────────────────────────────────────
    # mtr_index is 1-based → subtract 1 to get a 0-based component index.
    comp_idx = result.mtr_index - 1          # shape (N,), values 0..n_components-1
    comp_idx = np.clip(comp_idx, 0, n_components - 1)

    # Stack all samples into [n_components x N x 3], then index with comp_idx
    all_samples = np.stack(orient_samples, axis=0)           # (C, N, 3)
    voxel_idx   = np.arange(N)
    phi1_vec    = all_samples[comp_idx, voxel_idx, 0]
    phi_vec     = all_samples[comp_idx, voxel_idx, 1]
    phi2_vec    = all_samples[comp_idx, voxel_idx, 2]

    # ── 6. Write IPF map PNG ──────────────────────────────────────────────────
    ipf_path = os.path.join(args.output, "sim_IPF_map.png")
    print(f"Writing IPF map: {ipf_path}")
    mapper = mtrsim.IPFMapper(mtrsim.CrystalSystem.HCP)
    mapper.write_png(spatial_coords[:, :2],   # [N x 2] (x, y) only
                     phi1_vec, phi_vec, phi2_vec,
                     ipf_path)

    # ── 7. Write results CSV ──────────────────────────────────────────────────
    csv_path = os.path.join(args.output, "sim_results.csv")
    print(f"Writing CSV     : {csv_path}")
    write_csv(csv_path, spatial_coords, phi1_vec, phi_vec, phi2_vec, result.mtr_index)

    # ── 8. Quick numpy post-processing example ────────────────────────────────
    print()
    print("=== Post-processing (numpy) ===")

    # Mean misorientation spread (std dev of PHI) per component
    for comp_id in unique:
        mask    = result.mtr_index == comp_id
        phi_sub = phi_vec[mask]
        print(f"  Component {comp_id}: PHI mean={np.degrees(phi_sub.mean()):.2f}°  "
              f"std={np.degrees(phi_sub.std()):.2f}°  "
              f"n={mask.sum():,}")

    print()
    print("Done.")


if __name__ == "__main__":
    main()

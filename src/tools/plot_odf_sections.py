#!/usr/bin/env python3
"""
plot_odf_sections.py
Reads simulation_ODF.h5 and generates constant-phi2 section plots for each
component ODF. For HCP 6/mmm symmetry the fundamental zone is phi2 = 0..60°,
so we plot 13 sections at 5° intervals.

Usage:
    python plot_odf_sections.py <h5_path> <output_dir>
"""

import sys
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def load_odf_component(h5file, comp_idx):
    """Load a single component ODF from the HDF5 file."""
    grp = h5file[f"ODF_best/component_{comp_idx}"]
    odf_val = grp["ODFval"][:]
    phi1_bins = grp["phi1_bins"][:]
    phi_bins = grp["PHI_bins"][:]
    phi2_bins = grp["phi2_bins"][:]
    return odf_val, phi1_bins, phi_bins, phi2_bins


def reshape_odf(odf_val, n_phi1, n_PHI, n_phi2):
    """Reshape flat ODF array into 3D grid [phi1, PHI, phi2]."""
    # The ODF is stored as a flat array in the order:
    # phi1 varies fastest, then PHI, then phi2
    # Total bins = n_phi1 * n_PHI * n_phi2
    # But we need to figure out the actual storage order.
    # From ODFCalculator: ix = i2 * (72 * 36) + iPHI * 72 + i1
    # So: phi2 outer, PHI middle, phi1 inner
    odf_3d = odf_val.reshape(n_phi2, n_PHI, n_phi1)
    return odf_3d


def plot_component_odf(odf_3d, phi1_centres, phi_centres, phi2_centres,
                       comp_idx, output_path, global_vmax=None):
    """Plot constant-phi2 sections for one component ODF."""
    # HCP fundamental zone: phi2 = 0 to 60 degrees
    # phi2_centres are in radians; 60° = pi/3
    deg = np.degrees(phi2_centres)
    # Select sections in 0..60° range at 5° intervals
    target_angles = np.arange(0, 65, 5)  # 0, 5, 10, ..., 60
    section_indices = []
    for ta in target_angles:
        idx = np.argmin(np.abs(deg - ta))
        section_indices.append(idx)

    n_sections = len(section_indices)
    n_cols = 5
    n_rows = (n_sections + n_cols - 1) // n_cols

    # Use shared global max if provided, otherwise compute from this component
    if global_vmax is not None:
        global_max = global_vmax
    else:
        global_max = 0
        for si in section_indices:
            section = odf_3d[si, :, :]
            if section.max() > global_max:
                global_max = section.max()

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, n_rows * 3.2))
    fig.suptitle(f"ODF Component {comp_idx + 1} — Constant $\\varphi_2$ Sections (HCP 6/mmm)",
                 fontsize=14, fontweight='bold', y=0.98)

    phi1_deg = np.degrees(phi1_centres)
    phi_deg = np.degrees(phi_centres)

    norm = Normalize(vmin=0, vmax=global_max)

    for i, si in enumerate(section_indices):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]

        section = odf_3d[si, :, :]  # shape: [PHI, phi1]
        im = ax.pcolormesh(phi1_deg, phi_deg, section,
                           norm=norm, cmap='gray_r', shading='auto')
        ax.set_title(f"$\\varphi_2$ = {deg[si]:.0f}°", fontsize=10)
        ax.set_xlim(0, 360)
        ax.set_ylim(0, 180)
        ax.set_aspect('equal')
        ax.tick_params(labelsize=7)

        if col == 0:
            ax.set_ylabel("$\\Phi$ (°)", fontsize=9)
        if row == n_rows - 1:
            ax.set_xlabel("$\\varphi_1$ (°)", fontsize=9)

    # Hide unused axes
    for i in range(n_sections, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]
        ax.set_visible(False)

    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("ODF Intensity", fontsize=10)

    plt.subplots_adjust(left=0.06, right=0.90, top=0.92, bottom=0.08,
                        hspace=0.35, wspace=0.25)
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {output_path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: plot_odf_sections.py <h5_path> <output_dir>")
        sys.exit(1)

    h5_path = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)

    with h5py.File(h5_path, 'r') as f:
        num_components = int(f["ODF_best/num_components"][()])

        # First pass: load all ODFs and find global max across all components
        components = []
        global_max = 0.0
        for c in range(num_components):
            odf_val, phi1_bins, phi_bins, phi2_bins = load_odf_component(f, c)
            n_phi1 = len(phi1_bins) - 1
            n_PHI = len(phi_bins) - 1
            n_phi2 = len(phi2_bins) - 1
            phi1_centres = 0.5 * (phi1_bins[:-1] + phi1_bins[1:])
            phi_centres = 0.5 * (phi_bins[:-1] + phi_bins[1:])
            phi2_centres = 0.5 * (phi2_bins[:-1] + phi2_bins[1:])
            odf_3d = reshape_odf(odf_val, n_phi1, n_PHI, n_phi2)
            components.append((odf_3d, phi1_centres, phi_centres, phi2_centres))
            comp_max = odf_3d.max()
            if comp_max > global_max:
                global_max = comp_max

        print(f"Global ODF max across all components: {global_max:.6e}")

        # Second pass: plot with shared color scale
        for c, (odf_3d, phi1_centres, phi_centres, phi2_centres) in enumerate(components):
            print(f"Processing Component {c + 1}...")
            out_path = os.path.join(output_dir, f"odf_component_{c + 1}.png")
            plot_component_odf(odf_3d, phi1_centres, phi_centres, phi2_centres,
                               c, out_path, global_vmax=global_max)

    print("Done.")


if __name__ == "__main__":
    main()

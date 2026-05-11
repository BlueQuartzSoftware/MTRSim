#!/usr/bin/env python3
"""
plot_pole_figure.py — render a [0001] pole figure from sim_results.csv.

Usage:
    python tools/plot_pole_figure.py [results_csv] [output_png]

Defaults:
    results_csv : sim_results.csv  (in the current directory)
    output_png  : sim_pf.png       (in the same directory as results_csv)
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter


def load_euler_angles(csv_path):
    """Read phi1, PHI, phi2 columns from sim_results.csv using numpy."""
    print(f"Reading {csv_path} ...")
    # Identify column indices from header
    with open(csv_path) as f:
        header = f.readline().strip().split(",")
    col_phi1 = header.index("phi1")
    col_PHI  = header.index("PHI")
    col_phi2 = header.index("phi2")

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1,
                      usecols=(col_phi1, col_PHI, col_phi2))
    phi1 = data[:, 0]
    PHI  = data[:, 1]
    phi2 = data[:, 2]
    print(f"  Loaded {len(phi1):,} orientations.")
    return phi1, PHI, phi2


def compute_caxis(phi1, PHI):
    """
    [0001] c-axis in the sample frame.

    This is the last row of the Bunge passive rotation matrix G:
      h = [ sin(phi1)*sin(PHI),  -cos(phi1)*sin(PHI),  cos(PHI) ]

    Matches the convention in PoleFigure.cpp (convert_ODF_to_PF.m port).
    """
    hx = np.sin(phi1) * np.sin(PHI)
    hy = -np.cos(phi1) * np.sin(PHI)
    hz = np.cos(PHI)
    return hx, hy, hz


def stereographic_upper(hx, hy, hz):
    """
    Project to the upper hemisphere (hz >= 0) via south-pole stereographic:
      X = hx / (1 + hz),   Y = hy / (1 + hz)

    Points in the lower hemisphere are reflected through the origin first.
    """
    mask = hz < 0.0
    hx = hx.copy()
    hy = hy.copy()
    hz = hz.copy()
    hx[mask] = -hx[mask]
    hy[mask] = -hy[mask]
    hz[mask] = -hz[mask]

    denom = 1.0 + hz
    # Guard against exact south-pole singularity (hz == -1 before flip, now hz == 1)
    denom = np.where(denom < 1e-12, 1e-12, denom)
    X = hx / denom
    Y = hy / denom
    return X, Y


def build_density(X, Y, n_bins=200, sigma=2.5):
    """
    Bin (X, Y) into a 2-D histogram over [-1, 1]^2 and apply Gaussian
    smoothing.  Bins outside the unit circle are set to NaN.
    """
    H, xe, ye = np.histogram2d(X, Y, bins=n_bins, range=[[-1.0, 1.0], [-1.0, 1.0]])

    # Smooth
    H = gaussian_filter(H.astype(float), sigma=sigma)

    # Bin centres
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    Xg, Yg = np.meshgrid(xc, yc, indexing="ij")

    # Mask outside unit circle
    outside = (Xg ** 2 + Yg ** 2) > 1.0
    H[outside] = np.nan

    return Xg, Yg, H


def plot_pole_figure(Xg, Yg, H, out_path, title="[0001] Pole Figure"):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")

    vmin = np.nanmin(H)
    vmax = np.nanmax(H)
    n_levels = 14

    levels = np.linspace(vmin, vmax, n_levels)
    cf = ax.contourf(Xg, Yg, H, levels=levels, cmap="jet", extend="both")
    ax.contour(Xg, Yg, H, levels=levels, colors="k", linewidths=0.3, alpha=0.4)

    # Unit-circle boundary
    circle = Circle((0, 0), 1.0, fill=False, edgecolor="black", linewidth=1.5)
    ax.add_patch(circle)

    # Cardinal labels (RD/TD for a standard sample frame)
    offset = 1.08
    ax.text(0, offset, "RD", ha="center", va="bottom", fontsize=10)
    ax.text(0, -offset, "", ha="center", va="top", fontsize=10)
    ax.text(offset, 0, "TD", ha="left", va="center", fontsize=10)
    ax.text(-offset, 0, "", ha="right", va="center", fontsize=10)

    cbar = plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Counts (smoothed)", fontsize=9)

    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_title(title, fontsize=12)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "sim_results.csv"
    out_path = sys.argv[2] if len(sys.argv) > 2 else \
        os.path.join(os.path.dirname(os.path.abspath(csv_path)), "sim_pf.png")

    phi1, PHI, phi2 = load_euler_angles(csv_path)
    hx, hy, hz = compute_caxis(phi1, PHI)
    X, Y = stereographic_upper(hx, hy, hz)
    Xg, Yg, H = build_density(X, Y)
    plot_pole_figure(Xg, Yg, H, out_path)


if __name__ == "__main__":
    main()

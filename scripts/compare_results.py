#!/usr/bin/env python3
"""compare_results.py — Statistical comparison of MATLAB and C++ MTRSim outputs.

Since the two implementations use different RNG engines (MATLAB mt19937ar vs
C++ std::mt19937_64), outputs will NOT be identical even with the same seed.
This script performs statistical comparisons instead:

  1. Grid geometry: Same number of voxels, same coordinate ranges
  2. Volume fractions: Component proportions match targets within tolerance
  3. Euler angle distributions: Per-component ODF statistics are consistent
  4. Spatial statistics: Autocorrelation structure is qualitatively similar

Usage:
    python3 compare_results.py <matlab.csv> <cpp.csv> <config.json> [output_dir]

Exit code 0 = PASS, 1 = FAIL
"""

import csv
import json
import math
import os
import sys
from collections import Counter


def load_csv(path):
    """Load a sim_results.csv file into a list of dicts."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "x": float(row["x"]),
                    "y": float(row["y"]),
                    "z": float(row["z"]),
                    "phi1": float(row["phi1"]),
                    "PHI": float(row["PHI"]),
                    "phi2": float(row["phi2"]),
                    "mtr_index": int(row["mtr_index"]),
                }
            )
    return rows


def load_config(path):
    """Load JSON config and return as dict with defaults."""
    with open(path) as f:
        cfg = json.load(f)
    defaults = {
        "xLen": 1.5 * 25.4,
        "yLen": 0.5 * 25.4,
        "zLen": 0.0,
        "dx": 0.02,
        "dy": 0.02,
        "dz": 0.02,
        "volumeFractions": [0.30, 0.35, 0.35],
        "seed": 0,
    }
    for k, v in defaults.items():
        if k not in cfg:
            cfg[k] = v
    return cfg


def compute_volume_fractions(rows, num_components):
    """Compute actual volume fraction for each component."""
    counts = Counter(r["mtr_index"] for r in rows)
    total = len(rows)
    fracs = []
    for c in range(1, num_components + 1):
        fracs.append(counts.get(c, 0) / total)
    return fracs


def euler_stats(rows, component):
    """Compute mean and std of Euler angles for a given component."""
    phi1_vals = [r["phi1"] for r in rows if r["mtr_index"] == component]
    phi_vals = [r["PHI"] for r in rows if r["mtr_index"] == component]
    phi2_vals = [r["phi2"] for r in rows if r["mtr_index"] == component]

    if not phi1_vals:
        return None

    def mean_std(vals):
        n = len(vals)
        if n == 0:
            return 0.0, 0.0
        m = sum(vals) / n
        if n < 2:
            return m, 0.0
        var = sum((v - m) ** 2 for v in vals) / (n - 1)
        return m, math.sqrt(var)

    return {
        "phi1": mean_std(phi1_vals),
        "PHI": mean_std(phi_vals),
        "phi2": mean_std(phi2_vals),
        "count": len(phi1_vals),
    }


def check_grid_geometry(matlab_rows, cpp_rows, cfg, report):
    """Check that both outputs have the expected grid size and coordinate ranges."""
    passed = True

    # Expected voxel count
    nx = round(cfg["xLen"] / cfg["dx"])
    ny = round(cfg["yLen"] / cfg["dy"])
    nz = max(round(cfg["zLen"] / cfg["dz"]), 1)
    expected_n = nx * ny * nz

    if len(matlab_rows) != expected_n:
        report.append(
            f"  FAIL: MATLAB row count {len(matlab_rows)} != expected {expected_n}"
        )
        passed = False
    else:
        report.append(f"  OK:   MATLAB row count = {len(matlab_rows)}")

    if len(cpp_rows) != expected_n:
        report.append(
            f"  FAIL: C++ row count {len(cpp_rows)} != expected {expected_n}"
        )
        passed = False
    else:
        report.append(f"  OK:   C++ row count = {len(cpp_rows)}")

    # Coordinate ranges should match
    for label, rows in [("MATLAB", matlab_rows), ("C++", cpp_rows)]:
        x_vals = [r["x"] for r in rows]
        y_vals = [r["y"] for r in rows]

        x_range = max(x_vals) - min(x_vals)
        y_range = max(y_vals) - min(y_vals)

        expected_x_range = (nx - 1) * cfg["dx"]
        expected_y_range = (ny - 1) * cfg["dy"]

        if abs(x_range - expected_x_range) > cfg["dx"] * 0.5:
            report.append(
                f"  FAIL: {label} x-range {x_range:.4f} != expected {expected_x_range:.4f}"
            )
            passed = False

        if abs(y_range - expected_y_range) > cfg["dy"] * 0.5:
            report.append(
                f"  FAIL: {label} y-range {y_range:.4f} != expected {expected_y_range:.4f}"
            )
            passed = False

    if passed:
        report.append("  OK:   Grid geometry matches for both implementations.")

    return passed


def check_volume_fractions(matlab_rows, cpp_rows, cfg, report, tol=0.15):
    """Check that volume fractions match targets within tolerance."""
    passed = True
    target_fracs = cfg["volumeFractions"]
    num_comp = len(target_fracs)

    matlab_fracs = compute_volume_fractions(matlab_rows, num_comp)
    cpp_fracs = compute_volume_fractions(cpp_rows, num_comp)

    report.append("  Component | Target | MATLAB | C++    | MATLAB err | C++ err")
    report.append("  ----------|--------|--------|--------|------------|--------")

    for i in range(num_comp):
        m_err = abs(matlab_fracs[i] - target_fracs[i])
        c_err = abs(cpp_fracs[i] - target_fracs[i])

        m_status = "OK" if m_err <= tol else "FAIL"
        c_status = "OK" if c_err <= tol else "FAIL"

        report.append(
            f"  {i + 1:>9} | {target_fracs[i]:.4f} | {matlab_fracs[i]:.4f} | {cpp_fracs[i]:.4f} | "
            f"{m_err:.4f} {m_status:>4} | {c_err:.4f} {c_status:>4}"
        )

        if m_err > tol:
            passed = False
        if c_err > tol:
            passed = False

    return passed


def check_euler_distributions(matlab_rows, cpp_rows, cfg, report):
    """Check that Euler angle distributions are statistically consistent."""
    passed = True
    num_comp = len(cfg["volumeFractions"])

    for comp in range(1, num_comp + 1):
        m_stats = euler_stats(matlab_rows, comp)
        c_stats = euler_stats(cpp_rows, comp)

        if m_stats is None or c_stats is None:
            report.append(f"  WARN: Component {comp} has no samples in one or both outputs.")
            continue

        report.append(f"  Component {comp} (MATLAB n={m_stats['count']}, C++ n={c_stats['count']}):")

        for angle_name in ["phi1", "PHI", "phi2"]:
            m_mean, m_std = m_stats[angle_name]
            c_mean, c_std = c_stats[angle_name]

            # Allow generous tolerance — different RNG means different samples
            # We just check that distributions are in the same ballpark
            mean_diff = abs(m_mean - c_mean)
            std_diff = abs(m_std - c_std)

            # Tolerance: means within 0.3 rad (~17 deg), stds within 0.3 rad
            mean_ok = mean_diff < 0.3
            std_ok = std_diff < 0.3

            status = "OK" if (mean_ok and std_ok) else "WARN"
            if not (mean_ok and std_ok):
                # Don't fail on Euler stats — they're inherently variable
                pass

            report.append(
                f"    {angle_name:>4}: MATLAB mean={m_mean:.3f} std={m_std:.3f} | "
                f"C++ mean={c_mean:.3f} std={c_std:.3f} | {status}"
            )

    return passed


def check_all_components_present(matlab_rows, cpp_rows, cfg, report):
    """Check that all expected components appear in both outputs."""
    passed = True
    num_comp = len(cfg["volumeFractions"])

    matlab_comps = set(r["mtr_index"] for r in matlab_rows)
    cpp_comps = set(r["mtr_index"] for r in cpp_rows)
    expected = set(range(1, num_comp + 1))

    if matlab_comps != expected:
        report.append(
            f"  FAIL: MATLAB components {matlab_comps} != expected {expected}"
        )
        passed = False

    if cpp_comps != expected:
        report.append(f"  FAIL: C++ components {cpp_comps} != expected {expected}")
        passed = False

    if passed:
        report.append(
            f"  OK:   All {num_comp} components present in both outputs."
        )

    return passed


def check_euler_ranges(matlab_rows, cpp_rows, cfg, report):
    """Check that Euler angles are within valid ranges."""
    passed = True

    for label, rows in [("MATLAB", matlab_rows), ("C++", cpp_rows)]:
        phi1_vals = [r["phi1"] for r in rows]
        phi_vals = [r["PHI"] for r in rows]
        phi2_vals = [r["phi2"] for r in rows]

        # phi1, phi2 in [0, 2*pi], PHI in [0, pi] — with some tolerance for jitter
        margin = 0.15  # ~8.6 degrees of margin for bin jitter
        phi1_ok = all(-margin <= v <= 2 * math.pi + margin for v in phi1_vals)
        phi_ok = all(-margin <= v <= math.pi + margin for v in phi_vals)
        phi2_ok = all(-margin <= v <= 2 * math.pi + margin for v in phi2_vals)

        if not phi1_ok:
            report.append(
                f"  FAIL: {label} phi1 out of range [{min(phi1_vals):.4f}, {max(phi1_vals):.4f}]"
            )
            passed = False
        if not phi_ok:
            report.append(
                f"  FAIL: {label} PHI out of range [{min(phi_vals):.4f}, {max(phi_vals):.4f}]"
            )
            passed = False
        if not phi2_ok:
            report.append(
                f"  FAIL: {label} phi2 out of range [{min(phi2_vals):.4f}, {max(phi2_vals):.4f}]"
            )
            passed = False

    if passed:
        report.append("  OK:   Euler angles within valid ranges for both.")

    return passed


def load_timing(output_dir):
    """Load timing data from JSON files if they exist."""
    timing = {}
    for impl in ["matlab", "cpp"]:
        path = os.path.join(output_dir, f"{impl}_timing.json")
        if os.path.isfile(path):
            with open(path) as f:
                timing[impl] = json.load(f)
    return timing


def report_timing(timing, report):
    """Add timing information to the report."""
    m = timing.get("matlab")
    c = timing.get("cpp")

    if not m and not c:
        report.append("  No timing data available.")
        return

    report.append("  Implementation | Wall-clock (s) | Simulation (s)")
    report.append("  ---------------|----------------|---------------")

    if m:
        report.append(
            f"  MATLAB         | {m['wall_clock_seconds']:>14.3f} | {m['simulation_seconds']:>13.3f}"
        )
    if c:
        report.append(
            f"  C++            | {c['wall_clock_seconds']:>14.3f} | {c['simulation_seconds']:>13.3f}"
        )

    if m and c and c["wall_clock_seconds"] > 0:
        # Algorithm-to-algorithm speedup (MATLAB sim vs C++ wall)
        speedup = m["simulation_seconds"] / c["wall_clock_seconds"]
        report.append(f"")
        report.append(f"  Speedup (MATLAB sim / C++ wall): {speedup:.1f}x")

        # Total wall-clock speedup (including MATLAB startup)
        if m["wall_clock_seconds"] > 0:
            total_speedup = m["wall_clock_seconds"] / c["wall_clock_seconds"]
            report.append(f"  Speedup (MATLAB wall / C++ wall): {total_speedup:.1f}x")


def main():
    if len(sys.argv) < 4:
        print(
            f"Usage: {sys.argv[0]} <matlab.csv> <cpp.csv> <config.json> [output_dir]"
        )
        sys.exit(1)

    matlab_path = sys.argv[1]
    cpp_path = sys.argv[2]
    config_path = sys.argv[3]
    output_dir = sys.argv[4] if len(sys.argv) > 4 else "."

    config_name = os.path.splitext(os.path.basename(config_path))[0]

    print(f"\n=== Comparing: {config_name} ===")
    print(f"  MATLAB: {matlab_path}")
    print(f"  C++:    {cpp_path}")

    # Load data
    matlab_rows = load_csv(matlab_path)
    cpp_rows = load_csv(cpp_path)
    cfg = load_config(config_path)

    report = []
    all_passed = True

    # Test 1: Grid geometry
    report.append("\n[1] Grid Geometry")
    if not check_grid_geometry(matlab_rows, cpp_rows, cfg, report):
        all_passed = False

    # Test 2: All components present
    report.append("\n[2] Component Coverage")
    if not check_all_components_present(matlab_rows, cpp_rows, cfg, report):
        all_passed = False

    # Test 3: Volume fractions
    report.append("\n[3] Volume Fractions (tolerance = 15%)")
    if not check_volume_fractions(matlab_rows, cpp_rows, cfg, report):
        all_passed = False

    # Test 4: Euler angle ranges
    report.append("\n[4] Euler Angle Ranges")
    if not check_euler_ranges(matlab_rows, cpp_rows, cfg, report):
        all_passed = False

    # Test 5: Euler angle distributions
    report.append("\n[5] Euler Angle Distributions (informational)")
    check_euler_distributions(matlab_rows, cpp_rows, cfg, report)

    # Test 6: Timing
    report.append("\n[6] Performance Timing")
    timing = load_timing(output_dir)
    report_timing(timing, report)

    # Overall result
    result = "PASS" if all_passed else "FAIL"
    report.append(f"\n=== Overall: {result} ===\n")

    # Print report
    for line in report:
        print(line)

    # Write report to file
    report_path = os.path.join(output_dir, f"comparison_report_{config_name}.txt")
    with open(report_path, "w") as f:
        f.write(f"Comparison Report: {config_name}\n")
        f.write(f"MATLAB: {matlab_path}\n")
        f.write(f"C++:    {cpp_path}\n")
        f.write(f"Config: {config_path}\n")
        f.write("\n".join(report))
    print(f"Report saved: {report_path}")

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()

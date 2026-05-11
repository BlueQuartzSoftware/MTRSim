#!/usr/bin/env python3
"""
Convert the in-tree MATLAB v7.3 .mat ODF fixtures into the canonical single-
component HDF5 layout consumed by mtrsim::readODFMetadata / readODFComponents.

Inputs (under data/):
  - blank_ODF.mat   : flat at root. Contains /ODFval (1x186624), /ODFbins (3x186624),
                      /phi1_bins (73x1), /PHI_bins (37x1), /phi2_bins (73x1). All ODFval
                      values are zero.
  - uniformODF.mat  : nested under /uniformODF/ with the same dataset names. All ODFval
                      values equal 1 / 186624 (a normalised uniform ODF).

Outputs (under data/):
  - blank_ODF.h5    written with path prefix '/blank_ODF'   and num_components = 1.
  - uniform_ODF.h5  written with path prefix '/uniform_ODF' and num_components = 1.

The destination layout, for either file, is:

    <prefix>/num_components : int64 scalar
    <prefix>/component_0/ODFval    : float64 1-D (len = nphi1 * nPHI * nphi2)
    <prefix>/component_0/phi1_bins : float64 1-D (len = nphi1 + 1)
    <prefix>/component_0/PHI_bins  : float64 1-D (len = nPHI  + 1)
    <prefix>/component_0/phi2_bins : float64 1-D (len = nphi2 + 1)

Bin-edge arrays come straight from the source .mat (stored as Nx1 or 1xN),
flattened to 1-D. The source ODFbins (3xN coord triples) are discarded because
they are regeneratable from the edges.

This script is idempotent: re-running overwrites the outputs with byte-identical
contents (dataset order is fixed; gzip compression is disabled so we avoid any
compressor-version non-determinism; h5py writes the group header deterministically
given the same file-open sequence).
"""

from __future__ import annotations

import os
import sys
import numpy as np
import h5py


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data")


def _read_1d(dset: h5py.Dataset) -> np.ndarray:
    """Read an Nx1 or 1xN dataset as a contiguous 1-D float64 vector."""
    arr = np.array(dset[()], dtype=np.float64)
    return arr.flatten()


def _extract_source(src_path: str, group_name: str | None) -> dict:
    """Extract ODFval + bin edges from a .mat file.

    group_name is the top-level group under which the datasets live, or
    None if they sit at the root.
    """
    out: dict = {}
    with h5py.File(src_path, "r") as f:
        root = f[group_name] if group_name is not None else f
        odfval = _read_1d(root["ODFval"])
        phi1_bins = _read_1d(root["phi1_bins"])
        PHI_bins = _read_1d(root["PHI_bins"])
        phi2_bins = _read_1d(root["phi2_bins"])
    out["ODFval"] = odfval
    out["phi1_bins"] = phi1_bins
    out["PHI_bins"] = PHI_bins
    out["phi2_bins"] = phi2_bins
    return out


def _write_canonical(dst_path: str, prefix: str, data: dict) -> None:
    """Write a single-component ODF file at dst_path under `prefix`.

    Overwrites the destination. Order of dataset creation is fixed so reruns
    produce byte-identical files.
    """
    # Strip leading/trailing slashes from prefix for the h5py group path.
    clean = prefix.strip("/")
    nphi1 = int(data["phi1_bins"].size - 1)
    nPHI = int(data["PHI_bins"].size - 1)
    nphi2 = int(data["phi2_bins"].size - 1)

    expected_tuples = nphi1 * nPHI * nphi2
    if data["ODFval"].size != expected_tuples:
        raise RuntimeError(
            f"ODFval size {data['ODFval'].size} does not match "
            f"nphi1 * nPHI * nphi2 = {nphi1}*{nPHI}*{nphi2} = {expected_tuples}"
        )

    # Remove any existing file so runs are fully deterministic.
    if os.path.exists(dst_path):
        os.remove(dst_path)

    with h5py.File(dst_path, "w") as out:
        grp = out.create_group(clean)
        grp.create_dataset("num_components", data=np.int64(1))
        comp = grp.create_group("component_0")
        comp.create_dataset("ODFval", data=data["ODFval"].astype(np.float64, copy=False))
        comp.create_dataset("phi1_bins", data=data["phi1_bins"].astype(np.float64, copy=False))
        comp.create_dataset("PHI_bins", data=data["PHI_bins"].astype(np.float64, copy=False))
        comp.create_dataset("phi2_bins", data=data["phi2_bins"].astype(np.float64, copy=False))


def _summarise(label: str, src: str, dst: str, prefix: str, data: dict) -> None:
    nphi1 = int(data["phi1_bins"].size - 1)
    nPHI = int(data["PHI_bins"].size - 1)
    nphi2 = int(data["phi2_bins"].size - 1)
    odfval = data["ODFval"]
    nonzero = int(np.count_nonzero(odfval))
    unique_vals = np.unique(odfval)
    sample_val = float(odfval[0]) if odfval.size > 0 else 0.0
    print(f"[{label}]")
    print(f"  src     : {src}")
    print(f"  dst     : {dst}")
    print(f"  prefix  : {prefix}")
    print(f"  comps   : 1")
    print(f"  bins    : nphi1={nphi1} nPHI={nPHI} nphi2={nphi2} (tuples={nphi1*nPHI*nphi2})")
    print(f"  ODFval  : nonzero={nonzero}/{odfval.size}  first={sample_val:.6e}  unique_count={unique_vals.size}")


def main() -> int:
    blank_src = os.path.join(DATA_DIR, "blank_ODF.mat")
    uniform_src = os.path.join(DATA_DIR, "uniformODF.mat")
    blank_dst = os.path.join(DATA_DIR, "blank_ODF.h5")
    uniform_dst = os.path.join(DATA_DIR, "uniform_ODF.h5")

    for p in (blank_src, uniform_src):
        if not os.path.exists(p):
            print(f"ERROR: source fixture missing: {p}", file=sys.stderr)
            return 1

    # blank_ODF.mat: datasets at root
    blank = _extract_source(blank_src, group_name=None)
    # Sanity: per task spec, blank_ODF's ODFval must be all zeros.
    nonzero_blank = int(np.count_nonzero(blank["ODFval"]))
    if nonzero_blank != 0:
        print(
            f"WARNING: blank_ODF.mat ODFval has {nonzero_blank} non-zero values (expected 0). "
            f"Proceeding anyway.",
            file=sys.stderr,
        )
    _write_canonical(blank_dst, "/blank_ODF", blank)
    _summarise("blank_ODF", blank_src, blank_dst, "/blank_ODF", blank)

    # uniformODF.mat: datasets live under /uniformODF/
    uniform = _extract_source(uniform_src, group_name="uniformODF")
    # Sanity: every ODFval should equal 1 / 186624 ~ 5.36e-6
    total = float(uniform["ODFval"].size)
    expected_val = 1.0 / total if total > 0 else 0.0
    unique = np.unique(uniform["ODFval"])
    if unique.size != 1 or abs(unique[0] - expected_val) > 1e-12:
        print(
            f"WARNING: uniformODF.mat ODFval is not perfectly uniform. "
            f"unique_count={unique.size} first={uniform['ODFval'][0]:.6e} expected={expected_val:.6e}",
            file=sys.stderr,
        )
    _write_canonical(uniform_dst, "/uniform_ODF", uniform)
    _summarise("uniform_ODF", uniform_src, uniform_dst, "/uniform_ODF", uniform)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

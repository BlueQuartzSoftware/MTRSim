#!/usr/bin/env python3
"""
convert_simulation_ODF.py

Converts data/simulation_ODF.mat (MATLAB v5 format) to an HDF5 file
(data/simulation_ODF.h5) readable by h5py / HighFive in C++.

HDF5 layout produced:
  /ODF_best/
      num_components   -- scalar int64 dataset
      /component_0/
          ODFval       -- (186624,) float64
          phi1_bins    -- (73,)     float64
          PHI_bins     -- (37,)     float64
          phi2_bins    -- (73,)     float64
      /component_1/
          ...
      /component_2/
          ...

All arrays are stored as 1-D (squeezed) C-contiguous float64.
"""

import sys
from pathlib import Path
import numpy as np
import scipy.io
import h5py

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
INPUT_PATH = DATA_DIR / "simulation_ODF.mat"
OUTPUT_PATH = DATA_DIR / "simulation_ODF.h5"


def convert():
    if not INPUT_PATH.exists():
        print(f"ERROR: input file not found: {INPUT_PATH}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading  : {INPUT_PATH}")
    # squeeze_me=True collapses singleton dimensions and unwraps struct arrays
    # struct_as_record=False gives attribute-style access (elem.ODFval, etc.)
    mat = scipy.io.loadmat(str(INPUT_PATH), squeeze_me=True, struct_as_record=False)
    odf_best = mat["ODF_best"]  # now a 1-D ndarray of mat_struct objects
    num_components = odf_best.shape[0]
    fields = odf_best[0]._fieldnames
    print(f"  Found {num_components} ODF component(s), fields: {fields}")

    print(f"Writing  : {OUTPUT_PATH}")
    with h5py.File(OUTPUT_PATH, "w") as hf:
        grp = hf.create_group("ODF_best")
        grp.create_dataset("num_components", data=np.int64(num_components))

        for i in range(num_components):
            elem = odf_best[i]
            comp_grp = grp.create_group(f"component_{i}")

            for field in fields:
                arr = np.asarray(getattr(elem, field), dtype=np.float64).squeeze()
                comp_grp.create_dataset(field, data=arr)
                print(f"  component_{i}/{field}: shape={arr.shape}")

    # Quick verification round-trip
    print("\nVerification:")
    with h5py.File(OUTPUT_PATH, "r") as hf:
        nc = int(hf["ODF_best/num_components"][()])
        print(f"  num_components = {nc}")
        for i in range(nc):
            prefix = f"  component_{i}"
            for field in ("ODFval", "phi1_bins", "PHI_bins", "phi2_bins"):
                ds = hf[f"ODF_best/component_{i}/{field}"]
                print(f"  {prefix}/{field}: shape={ds.shape} dtype={ds.dtype}")

    print(f"\nDone. Output written to: {OUTPUT_PATH}")


if __name__ == "__main__":
    convert()

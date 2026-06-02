#!/usr/bin/env python3
"""Generate DREAM3D-NX .d3dpipeline files from the LibMTRSim configs/*.json files.

Each generated pipeline reproduces one MTRSim configuration as a three-step
DREAM3D-NX pipeline:

    ReadMTRSimODFFilter  ->  MTRSimFilter  ->  WriteDREAM3DFilter

The MTRSim simulation parameters (volume fractions, theta list, physical size /
spacing, seed) are taken directly from the JSON config so the NX pipeline matches
the standalone LibMTRSim run. The number of ODF component arrays is derived from
the length of ``volumeFractions``.

Units note: the JSON configs are millimeter-consistent (size, spacing, and theta
all in mm). The MTRSimFilter "Physical Size/Spacing (microns)" labels are nominal
— the simulation only uses the dimensionless lag/theta ratio — so the mm values
are written through unchanged.

Usage:
    python3 scripts/generate_nx_pipelines.py [--configs-dir DIR] [--out-dir DIR]
"""

import argparse
import json
from pathlib import Path

# Filter identity constants (kept in sync with the filter headers).
READ_ODF_UUID = "2b1a4841-65d7-4315-9fe3-d66c88e5755c"
MTRSIM_UUID = "f7f7a330-4bff-4a42-a573-09117a89a0a0"
WRITE_DREAM3D_UUID = "b3a95784-2ced-41ec-8d3d-0242ac130003"

ODF_GEOM = "ODF"
ODF_CELL_AM = "Cell Data"


def _val(value, version=1):
    return {"value": value, "version": version}


def build_pipeline(name, cfg, repo_root):
    """Build the pipeline dict for one config."""
    vf = cfg["volumeFractions"]
    num_components = len(vf)
    component_paths = [f"{ODF_GEOM}/{ODF_CELL_AM}/component_{i}" for i in range(num_components)]

    odf_file = (repo_root / "data" / "simulation_ODF.h5").as_posix()
    out_file = (repo_root / "output" / f"MTRSim_{name}.dream3d").as_posix()

    read_step = {
        "args": {
            "input_file": _val(odf_file),
            "hdf5_path_prefix": _val("/ODF_best"),
            "output_image_geometry_path": _val(ODF_GEOM),
            "cell_attribute_matrix_name": _val(ODF_CELL_AM),
            "parameters_version": 1,
        },
        "comments": f"Load the {num_components}-component MTRSim ODF into ImageGeom '{ODF_GEOM}'.",
        "filter": {"name": "nx::core::ReadMTRSimODFFilter", "uuid": READ_ODF_UUID},
        "isDisabled": False,
    }

    mtrsim_step = {
        "args": {
            "input_odf_geometry_path": _val(ODF_GEOM),
            "odf_component_arrays": _val(component_paths),
            "volume_fractions": _val([list(vf)]),
            "theta_list": _val([list(row) for row in cfg["thetaList"]]),
            "physical_size": _val([cfg["xLen"], cfg["yLen"], cfg["zLen"]]),
            "physical_spacing": _val([cfg["dx"], cfg["dy"], cfg["dz"]]),
            "use_seed": _val(True),
            "seed_value": _val(int(cfg.get("seed", 42))),
            "seed_array_name": _val("MTRSim SeedValue"),
            "generate_polar_coloring": _val(True),
            "output_geometry_path": _val("MTR Microstructure"),
            "cell_attribute_matrix_name": _val("Cell Data"),
            "mtr_ids_array_name": _val("MTRIds"),
            "eulers_array_name": _val("Eulers"),
            "polar_colors_array_name": _val("Polar Colors"),
            "parameters_version": 1,
        },
        "comments": cfg.get("_comment", f"Recreates configs/{name}.json."),
        "filter": {"name": "nx::core::MTRSimFilter", "uuid": MTRSIM_UUID},
        "isDisabled": False,
    }

    write_step = {
        "args": {
            "export_file_path": _val(out_file),
            "write_xdmf_file": _val(False),
            "use_compression": _val(True),
            "compression_level": _val(5),
            "parameters_version": 1,
        },
        "comments": "Save the synthetic microstructure. Disable to only view in the GUI.",
        "filter": {"name": "nx::core::WriteDREAM3DFilter", "uuid": WRITE_DREAM3D_UUID},
        "isDisabled": False,
    }

    return {
        "isDisabled": False,
        "name": f"MTRSim_{name}.d3dpipeline",
        "pinnedParams": [],
        "pipeline": [read_step, mtrsim_step, write_step],
        "version": 1,
        "workflowParams": [],
    }


def main():
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configs-dir", default=str(repo_root / "configs"))
    parser.add_argument("--out-dir", default=str(repo_root / "pipelines"))
    args = parser.parse_args()

    configs_dir = Path(args.configs_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for cfg_path in sorted(configs_dir.glob("*.json")):
        with open(cfg_path) as f:
            cfg = json.load(f)
        if "volumeFractions" not in cfg:
            continue
        name = cfg_path.stem
        pipeline = build_pipeline(name, cfg, repo_root)
        out_path = out_dir / f"MTRSim_{name}.d3dpipeline"
        with open(out_path, "w") as f:
            json.dump(pipeline, f, indent=2)
            f.write("\n")
        count += 1
        print(f"wrote {out_path.relative_to(repo_root)}")
    print(f"\nGenerated {count} pipeline(s).")


if __name__ == "__main__":
    main()

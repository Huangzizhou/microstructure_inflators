#!/usr/bin/env python

import argparse
import json
import os.path

def output_config_file(wire_list_file, guide_mesh_file, thickness_type,
        dof_type, output_file, spread_const, abs_geometry_correction):
    root_dir = os.path.dirname(output_file);
    wire_list_file = os.path.relpath(wire_list_file, root_dir);
    guide_mesh_file = os.path.relpath(guide_mesh_file, root_dir);
    config = {
            "wire_list_file": wire_list_file,
            "geometry_spread": spread_const,
            "abs_geometry_correction": abs_geometry_correction,
            "guide_mesh": guide_mesh_file,
            "subdiv": 0,
            "dof_type": dof_type,
            "thickness_type": thickness_type
            };
    with open(output_file, 'w') as fout:
        json.dump(config, fout, indent=4);

def parse_args():
    parser = argparse.ArgumentParser(description="Generate config file");
    parser.add_argument("--thickness-type", choices=("vertex", "edge"),
            default="edge");
    parser.add_argument("--dof-type", choices=("isotropic", "orthotropic"),
            default="isotropic");
    parser.add_argument("--spread", help="spread constant", default=0,
            type=float);
    parser.add_argument("--abs-geometry-correction",
            help="absolute geometry correction", nargs=3, type=float,
            default=[0, 0, 0]);
    parser.add_argument("wire_list_file", help="wire list file");
    parser.add_argument("guide_mesh_file", help="guide mesh file");
    parser.add_argument("output_file", help="output config file");
    return parser.parse_args();

def main():
    args = parse_args();
    output_config_file(
            args.wire_list_file,
            args.guide_mesh_file,
            args.thickness_type,
            args.dof_type,
            args.output_file,
            args.spread,
            args.abs_geometry_correction);

if __name__ == "__main__":
    main();

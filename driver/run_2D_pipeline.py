#!/usr/bin/env python

""" Run 2D pipeline """

import argparse
import hashlib
import json
import os
import os.path
import re
from subprocess import check_call
import sys

PROJECT_DIR = os.environ["MICROSTRUCTURES_PATH"];

def get_stamp(name):
    return hashlib.md5(name).hexdigest();

def generate_config(wire_file, quad_mesh_file, modifier_file):
    config = {
            "wire_network": wire_file,
            "thickness": 0.5,
            "modifier_file": modifier_file,
            "guide_mesh": quad_mesh_file,
            }
    return config;

def parse_args():
    parser = argparse.ArgumentParser(description="Run all steps of 2D pipeline");
    parser.add_argument("--output", "-o", help="output microstructure", required=True);
    parser.add_argument("--pattern", help="pattern name", default="box_2D");
    parser.add_argument("--cell-size", help="size of each cell", type=float,
            default=0.5);
    parser.add_argument("--metric", help="material metric",
            choices=("compliance", "elasticity"));
    parser.add_argument("quad_mesh", help="quad mesh defining cell layout");
    parser.add_argument("material_mesh",
            help="material mesh containt material property attributes");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    spliter = re.compile(r'(?:[^\s,"]|"(?:\\.|[^"])*")+');

    TMP_DIR = "/tmp"

    basename, ext = os.path.splitext(args.quad_mesh);
    path, name = os.path.split(basename);
    stamp = get_stamp(name);

    # Copy quad mesh to a tmp dir.
    tmp_material_mesh = os.path.join(TMP_DIR, stamp + "_quad.msh");
    command = "cp {} {}".format(args.quad_mesh, tmp_material_mesh);
    check_call(command.split());

    # Rename "Final E" and "Final nu" to "young" and "poisson"
    exe_name = os.path.join(PROJECT_DIR, "tools/msh_tools/rename_element_attribute.py");
    command = [exe_name, args.material_mesh, "Final E", "young"];
    check_call(command);
    command = [exe_name, args.material_mesh, "Final nu", "poisson"];
    check_call(command);

    # Copy material attributes from material mesh to quad mesh.
    exe_name = os.path.join(PROJECT_DIR, "tools/msh_tools/move_element_attribute.py");
    command = "{} --attribute young --attribute poisson {} {}".format(exe_name, args.material_mesh,
            tmp_material_mesh)
    check_call(command.split());

    # Convert material parameter to pattern parameter
    tmp_pattern_mesh = os.path.join(TMP_DIR, stamp + "_pattern.msh");
    index_dir = os.path.join(PROJECT_DIR,
            "lookup_table/configs/2D/{}/{}mm_cell/mixed/index"\
                    .format(args.pattern, int(args.cell_size)));
    assert(os.path.isdir(index_dir));
    exe_name = os.path.join(PROJECT_DIR, "parameter_lookup/lookup.py");
    command = "{} --metric {} --index-dir {} {} {}".format(exe_name,
            args.metric, index_dir, tmp_material_mesh, tmp_pattern_mesh);
    check_call(command.split());

    # Tile microstructure according to pattern parameter
    wire_file = "patterns/2D/{}.wire".format(args.pattern);
    exe_name = os.path.join(PROJECT_DIR, "wire_generator/tile.py");
    modifier_file = os.path.join(PROJECT_DIR,
            "{}/lookup.modifier".format(index_dir));
    modifier_file = os.path.abspath(modifier_file);
    config = generate_config(wire_file,
            os.path.abspath(tmp_pattern_mesh), modifier_file);
    tmp_config_file = os.path.join(TMP_DIR, stamp + ".config");
    with open(tmp_config_file, 'w') as fout:
        json.dump(config, fout);
    command = "{} -o {} {}".format(exe_name, args.output, tmp_config_file);
    check_call(command.split());

    # Move pattern mesh over.
    basename, ext = os.path.splitext(args.output);
    pattern_mesh_filename = basename + "_reference.msh";
    os.rename(tmp_pattern_mesh, pattern_mesh_filename);

    os.remove(tmp_material_mesh);
    #os.remove(tmp_pattern_mesh);
    os.remove(tmp_config_file);

if __name__ == "__main__":
    main();



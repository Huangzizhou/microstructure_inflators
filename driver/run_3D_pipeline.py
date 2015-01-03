#!/usr/bin/env python

""" Run 3D pipeline """

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

def generate_config(wire_file, guide_mesh_file, modifier_file):
    config = {
            "trim": True,
            "wire_network": wire_file,
            "subdiv": 0,
            "thickness": 0.5,
            "modifier_file": modifier_file,
            "guide_mesh": guide_mesh_file,
            }
    return config;

def parse_args():
    parser = argparse.ArgumentParser(description="Run all steps of 3D pipeline");
    parser.add_argument("--output", "-o", help="output microstructure", required=True);
    parser.add_argument("--pattern", help="pattern name",
            default="truncated_octahedron");
    parser.add_argument("--cell-size", help="size of each cell", type=str,
            default="5");
    parser.add_argument("--sweep-type", help="sweep type",
            choices=("thickness_only", "offset_only", "mixed", "isotropic"),
            default="thickness_only");
    parser.add_argument("--metric", help="material metric",
            choices=("compliance", "elasticity"));
    parser.add_argument("--rehomogenize", action="store_true",
            help="run homogenization on interpolated pattern parameters");
    parser.add_argument("--material", default=None,
            help="material filed used");
    parser.add_argument("guide_mesh", help="guide mesh defining cell layout");
    parser.add_argument("material_mesh",
            help="material mesh containt material property attributes");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    spliter = re.compile(r'(?:[^\s,"]|"(?:\\.|[^"])*")+');

    TMP_DIR = "/tmp"

    basename, ext = os.path.splitext(args.guide_mesh);
    path, name = os.path.split(basename);
    stamp = get_stamp(name);

    # Copy hex mesh to a tmp dir.
    tmp_material_mesh = os.path.join(TMP_DIR, stamp + "_guide.msh");
    command = "cp {} {}".format(args.guide_mesh, tmp_material_mesh);
    check_call(command.split());

    # Rename "Final E" and "Final nu" to "young" and "poisson"
    exe_name = os.path.join(PROJECT_DIR, "tools/msh_tools/rename_element_attribute.py");
    command = [exe_name, args.material_mesh, "Final E", "young"];
    check_call(command);
    command = [exe_name, args.material_mesh, "Final nu", "poisson"];
    check_call(command);

    # Copy material attributes from material mesh to guide mesh.
    exe_name = os.path.join(PROJECT_DIR, "tools/msh_tools/move_element_attribute.py");
    command = "{} --attribute young --attribute poisson {} {}".format(exe_name, args.material_mesh,
            tmp_material_mesh)
    check_call(command.split());

    # Convert material parameter to pattern parameter
    tmp_pattern_mesh = os.path.join(TMP_DIR, stamp + "_pattern.msh");
    index_dir = os.path.join(PROJECT_DIR,
            "lookup_table/configs/3D/{}/{}mm_cell/{}/index"\
                    .format(args.pattern, args.cell_size, args.sweep_type));
    assert(os.path.isdir(index_dir));
    exe_name = os.path.join(PROJECT_DIR, "parameter_lookup/lookup.py");
    if args.rehomogenize:
        if args.material is None or not os.path.exists(args.material):
            raise RuntimeError("Material file ({}) is invalid".format(
                args.material));
        command = "{} --rehomogenize --material {} --metric {} --index-dir {} {} {}".format(
                exe_name, args.material, args.metric, index_dir,
                tmp_material_mesh, tmp_pattern_mesh);
    else:
        command = "{} --metric {} --index-dir {} {} {}".format(exe_name,
                args.metric, index_dir, tmp_material_mesh, tmp_pattern_mesh);
    check_call(command.split());

    # Tile microstructure according to pattern parameter
    wire_file = "patterns/3D/{}.wire".format(args.pattern);
    exe_name = os.path.join(PROJECT_DIR, "wire_generator/tile.py");
    config_file = os.path.join(TMP_DIR, stamp + "_pattern.config");
    #modifier_file = os.path.join(PROJECT_DIR,
    #        "{}/lookup.modifier".format(index_dir));
    #modifier_file = os.path.abspath(modifier_file);
    #config = generate_config(wire_file,
    #        os.path.abspath(tmp_pattern_mesh), modifier_file);
    #tmp_config_file = os.path.join(TMP_DIR, stamp + ".config");
    #with open(tmp_config_file, 'w') as fout:
    #    json.dump(config, fout);
    command = "{} -o {} {}".format(exe_name, args.output, config_file);
    check_call(command.split());

    # Move pattern mesh over.
    basename, ext = os.path.splitext(args.output);
    pattern_mesh_filename = basename + "_reference.msh";
    os.rename(tmp_pattern_mesh, pattern_mesh_filename);
    #config_file = basename + ".config";
    #os.rename(tmp_config_file, config_file);

    os.remove(tmp_material_mesh);
    #os.remove(tmp_pattern_mesh);
    #os.remove(tmp_config_file);

if __name__ == "__main__":
    main();



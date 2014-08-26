#!/usr/bin/env python

import argparse
import json
import hashlib
from subprocess import check_call, check_output
import os
import os.path
import re
import sys

def generate_stamp(config_file):
    m = hashlib.md5();
    m.update(config_file);
    return m.hexdigest();

def run_tile(config_file, obj_file):
    exe_dir = sys.path[0];
    exe_name = os.path.join(exe_dir, "tile.py");
    cmd = "{} -o {} {}".format(exe_name, obj_file, config_file);
    print(cmd);
    check_call(cmd.split());

def run_material_fit(input_file, material_file, output_file):
    meshfem_path = os.environ.get("MESHFEM_PATH");
    exe_name = os.path.join(meshfem_path, "PeriodicHomogenization_cli");
    cmd = "{} --outputMSH {} --material {} {}".format(exe_name,
            output_file, material_file, input_file);
    print(cmd);
    result = check_output(cmd.split());

    young, shear, poisson = parse_result(result);
    mat_config = {
            "youngs_modulus": young,
            "shear_modulus": shear,
            "poisson_ratio": poisson };

    basename, ext = os.path.splitext(output_file);
    json_file = "{}_param.json".format(basename);
    with open(json_file, 'w') as fout:
        json.dump(mat_config, fout, indent=4);
    print("material property written in {}".format(json_file));

def parse_young(result):
    young_pattern = "Approximate Young moduli:\s*(\d+\.?\d+)\s*(\d+\.?\d+)";
    young_matcher = re.compile(young_pattern, re.M);
    young_result = young_matcher.search(result);
    assert(young_result is not None);

    young_x = float(young_result.group(1));
    young_y = float(young_result.group(2));
    return [young_x, young_y];

def parse_shear(result):
    shear_pattern = "Approximate shear modulus:\s*(\d+\.?\d+)";
    shear_matcher = re.compile(shear_pattern, re.M);
    shear_result = shear_matcher.search(result);
    assert(shear_result is not None);

    shear_xy = float(shear_result.group(1));
    return [shear_xy];

def parse_poisson(result):
    poisson_1_pattern = "v_yx, v_xy:\s*(\d+\.?\d+)\s*(\d+\.?\d+)";
    poisson_1_matcher = re.compile(poisson_1_pattern, re.M);
    poisson_1_result = poisson_1_matcher.search(result);
    assert(poisson_1_result is not None);

    v_yx = float(poisson_1_result.group(1));
    v_xy = float(poisson_1_result.group(2));

    return [v_xy, v_yx];

def parse_result(result):
    young = parse_young(result);
    shear = parse_shear(result);
    poisson = parse_poisson(result);

    return young, shear, poisson;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Tile pattern, tetgen it, then fit orthotropic material");
    parser.add_argument("--material", help="material file");
    parser.add_argument("config_file", help="configuration file");
    parser.add_argument("msh_file", help="output msh file");
    args = parser.parse_args();
    return args

def main():
    args = parse_args();

    stamp = generate_stamp(args.msh_file);
    tmp_dir = "/tmp"
    tmp_msh = os.path.join(tmp_dir, stamp+".msh");

    run_tile(args.config_file, tmp_msh);
    run_material_fit(tmp_msh, args.material, args.msh_file);

    os.remove(tmp_msh);

if __name__ == "__main__":
    main();


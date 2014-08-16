#!/usr/bin/env python

import argparse
import json
import hashlib
from subprocess import check_call, check_output
import os
import os.path
import re

def generate_stamp(config_file):
    m = hashlib.md5();
    m.update(config_file);
    return m.hexdigest();

def run_tile(config_file, obj_file):
    cmd = "./tile.py -o {} {}".format(obj_file, config_file);
    print(cmd);
    check_call(cmd.split());

def run_tetgen(obj_file, msh_file):
    cmd = "tetgen.py --cmd --flags=\"qpQa0.5\" {} {}".format(obj_file, msh_file);
    print(cmd);
    check_call(cmd.split());

def run_material_fit(input_file, material_file, output_file):
    meshfem_path = os.environ.get("MESHFEM_PATH");
    exe_name = os.path.join(meshfem_path, "PeriodicHomogenization_cli");
    cmd = "{} --material {} {}".format(exe_name,
            material_file, input_file);
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
    young_pattern = "Approximate Young moduli:\s*(\d+\.?\d+)\s*(\d+\.?\d+)\s*(\d+\.?\d+)";
    young_matcher = re.compile(young_pattern, re.M);
    young_result = young_matcher.search(result);
    assert(young_result is not None);

    young_x = float(young_result.group(1));
    young_y = float(young_result.group(2));
    young_z = float(young_result.group(3));
    return [young_x, young_y, young_z];

def parse_shear(result):
    shear_pattern = "Approximate shear moduli:\s*(\d+\.?\d+)\s*(\d+\.?\d+)\s*(\d+\.?\d+)";
    shear_matcher = re.compile(shear_pattern, re.M);
    shear_result = shear_matcher.search(result);
    assert(shear_result is not None);

    shear_yz = float(shear_result.group(1));
    shear_zx = float(shear_result.group(2));
    shear_xy = float(shear_result.group(3));
    return [shear_yz, shear_zx, shear_xy];

def parse_poisson(result):
    poisson_1_pattern = "v_yx, v_zx, v_zy:\s*(\d+\.?\d+)\s*(\d+\.?\d+)\s*(\d+\.?\d+)";
    poisson_1_matcher = re.compile(poisson_1_pattern, re.M);
    poisson_1_result = poisson_1_matcher.search(result);
    assert(poisson_1_result is not None);

    poisson_2_pattern = "v_xy, v_xz, v_yz:\s*(\d+\.?\d+)\s*(\d+\.?\d+)\s*(\d+\.?\d+)";
    poisson_2_matcher = re.compile(poisson_2_pattern, re.M);
    poisson_2_result = poisson_2_matcher.search(result);
    assert(poisson_2_result is not None);

    v_yx = float(poisson_1_result.group(1));
    v_zx = float(poisson_1_result.group(2));
    v_zy = float(poisson_1_result.group(3));

    v_xy = float(poisson_2_result.group(1));
    v_xz = float(poisson_2_result.group(2));
    v_yz = float(poisson_2_result.group(3));

    return [v_yz, v_zy, v_zx, v_xz, v_xy, v_yx];

def parse_result(result):
    young = parse_young(result);
    shear = parse_shear(result);
    poisson = parse_poisson(result);

    return young, shear, poisson;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Tile pattern, tetgen it, then fit orthotropic material");
    parser.add_argument("--material", help="material file", required=True);
    parser.add_argument("config_file", help="configuration file");
    parser.add_argument("msh_file", help="output msh file");
    args = parser.parse_args();
    return args

def main():
    args = parse_args();

    stamp = generate_stamp(args.msh_file);
    tmp_dir = "/tmp"
    tmp_obj = os.path.join(tmp_dir, stamp+".obj");
    #tmp_msh = os.path.join(tmp_dir, stamp+".msh");

    run_tile(args.config_file, tmp_obj);
    run_tetgen(tmp_obj, args.msh_file);
    run_material_fit(args.msh_file, args.material, args.msh_file);

    os.remove(tmp_obj);
    #os.remove(tmp_msh);

if __name__ == "__main__":
    main();


#!/usr/bin/env python

import argparse
import hashlib
from subprocess import check_call
import os
import os.path

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

def run_orthotropic_fit(input_file, material_file, output_file):
    orthotropic_fit_path = os.environ.get("ORTHOTROPIC_FIT_PATH");
    exe_name = os.path.join(orthotropic_fit_path, "fit_orthotropic_material.py");
    cmd = "{} --material {} {} {}".format(exe_name, material_file, input_file,
            output_file);
    print(cmd);
    check_call(cmd.split());

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
    tmp_obj = os.path.join(tmp_dir, stamp+".obj");
    tmp_msh = os.path.join(tmp_dir, stamp+".msh");

    run_tile(args.config_file, tmp_obj);
    run_tetgen(tmp_obj, tmp_msh);
    run_orthotropic_fit(tmp_msh, args.material, args.msh_file);

    os.remove(tmp_obj);
    os.remove(tmp_msh);

if __name__ == "__main__":
    main();


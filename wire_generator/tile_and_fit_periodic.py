#!/usr/bin/env python

import argparse
import json
import numpy as np
from numpy.linalg import eig
import os
import os.path
import re
from subprocess import check_call, check_output

import core.PyMeshSetting
import PyMesh
from tile import tile, parse_config_file
from tile_and_tetgen import generate_tetgen_flags, run_tetgen
from utils.mesh_io import save_mesh_raw
from utils.timethis import timethis
from utils.print_colors import bcolors

@timethis
def run_material_fit(input_file, material_file, output_file):
    meshfem_path = os.environ.get("MESHFEM_PATH");
    exe_name = os.path.join(meshfem_path, "PeriodicHomogenization_cli");
    cmd = "{} --material {} {}".format(exe_name,
            material_file, input_file);
    print(cmd);
    result = check_output(cmd.split());

    young, shear, poisson, elasticity = parse_result(result);
    modes = compute_modes(elasticity);
    mat_config = {
            "youngs_modulus": young,
            "shear_modulus": shear,
            "poisson_ratio": poisson,
            "elasticity_tensor": elasticity,
            "elasticity_modes": modes};

    basename, ext = os.path.splitext(output_file);
    json_file = "{}_param.json".format(basename);
    with open(json_file, 'w') as fout:
        json.dump(mat_config, fout, indent=4);
    print("material property written in {}".format(json_file));

def parse_field(input_text, header, num_entries):
    pattern = "^{}(.*)$".format(header);
    matcher = re.compile(pattern, re.M);
    match_result = matcher.search(input_text);
    assert(match_result is not None);
    result = map(float, match_result.group(1).split());
    assert(len(result) == num_entries);
    return result;

def parse_multiline_field(input_text, header, tailer, num_entries):
    pattern = "{}(.*){}".format(header, tailer);
    matcher = re.compile(pattern, re.M | re.S);
    match_result = matcher.search(input_text);
    assert(match_result is not None);
    result = map(float, match_result.group(1).split());
    assert(len(result) == num_entries);
    return result;

def parse_young(result):
    header = "Approximate Young moduli:";
    young_x, young_y, young_z = parse_field(result, header, 3);
    return [young_x, young_y, young_z];

def parse_shear(result):
    header = "Approximate shear moduli:"
    shear_yz, shear_zx, shear_xy = parse_field(result, header, 3);
    return [shear_yz, shear_zx, shear_xy];

def parse_poisson(result):
    header1 = "v_yx, v_zx, v_zy:";
    header2 = "v_xy, v_xz, v_yz:";
    v_yx, v_zx, v_zy = parse_field(result, header1, 3);
    v_xy, v_xz, v_yz = parse_field(result, header2, 3);
    return [v_yz, v_zy, v_zx, v_xz, v_xy, v_yx];

def parse_elasticity_tensor(result):
    header = "Homogenized elasticity tensor:";
    tailer = "\n\n";
    elasticity = parse_multiline_field(result, header, tailer, 36);
    return np.reshape(elasticity, (6, 6)).tolist();

def parse_result(result):
    young = parse_young(result);
    shear = parse_shear(result);
    poisson = parse_poisson(result);
    elasticity = parse_elasticity_tensor(result);

    return young, shear, poisson, elasticity;

def compute_modes(elasticity_tensor):
    eig_val, eig_vec = eig(elasticity_tensor);
    eig_val = sorted(eig_val, reverse=True);
    if (eig_val[0] > 100 * eig_val[1]):
        print(bcolors.OKGREEN + "This pattern might be pentamode!!!" +
                bcolors.ENDC);
    return eig_val;

def save_timing(msh_file):
    basename, ext = os.path.splitext(msh_file);
    timing_file = "{}.timing".format(basename);
    timer = timethis(None);
    with open(timing_file, 'w') as fout:
        json.dump(timer.hist, fout, indent=4);

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

    config = parse_config_file(args.config_file);
    mesh = tile(config);
    vertices, faces, voxels = run_tetgen(config, mesh);
    save_mesh_raw(args.msh_file, vertices, faces, voxels);
    run_material_fit(args.msh_file, args.material, args.msh_file);

    timethis.summarize();
    save_timing(args.msh_file);

if __name__ == "__main__":
    main();


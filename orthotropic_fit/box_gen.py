#!/usr/bin/env python

import argparse
import numpy as np

import LinearElasticitySettings
from mesh_io import save_mesh, load_mesh

from BoxMeshGenerator import generate_box_mesh

def generate_box(dim, side_length, num_samples, output_name, keep_symmetry):
    box_min = -0.5 * np.ones(dim) * side_length;
    box_max =  0.5 * np.ones(dim) * side_length;
    mesh = generate_box_mesh(box_min, box_max, num_samples, keep_symmetry);
    index = np.arange(num_samples**3, dtype=int)
    if dim == 2:
        index = np.repeat(index, 2);
    elif keep_symmetry:
        index = np.repeat(index, 24);
    else:
        index = np.repeat(index, 6);
    mesh.add_attribute("hex_index");
    mesh.set_attribute("hex_index", index);
    save_mesh(output_name, mesh, "hex_index");

def parse_args():
    parser = argparse.ArgumentParser(description="Generate box meshes");
    parser.add_argument("-s", "--symmetric", help="symmetric tet connectivity",
            action="store_true");
    parser.add_argument("--dim", help="mesh dimention", choices=[2, 3],
            default=3, type=int);
    parser.add_argument("--size", help="box size", type=float, default=10);
    parser.add_argument("--num-samples",
            help="number of samples along each dimention", type=int, default=2);
    parser.add_argument("output", help="output_mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    generate_box(args.dim, args.size, args.num_samples, args.output,
            args.symmetric);

if __name__ == "__main__":
    main();

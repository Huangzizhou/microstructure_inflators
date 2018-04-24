#!/usr/bin/env python

import argparse
import numpy as np

import LinearElasticitySettings
from mesh_io import save_mesh, load_mesh

from BoxMeshGenerator import generate_box_mesh

def generate_box(side_lengths, num_samples, output_name, keep_symmetry,
        subdiv_order):
    box_min = -0.5 * side_lengths;
    box_max =  0.5 * side_lengths;
    mesh, material_index = generate_box_mesh(box_min, box_max, num_samples, keep_symmetry,
            subdiv_order);
    index_name = "cell_index"
    mesh.add_attribute(index_name);
    mesh.set_attribute(index_name, material_index);
    save_mesh(output_name, mesh, index_name);

def parse_args():
    parser = argparse.ArgumentParser(description="Generate box meshes");
    parser.add_argument("-s", "--symmetric", help="symmetric tet connectivity",
            action="store_true");
    parser.add_argument("--dim", help="mesh dimention", choices=[2, 3],
            default=3, type=int);
    parser.add_argument("--size", help="box size", type=float, default=10);
    parser.add_argument("-X", help="X size", type=float, default=None);
    parser.add_argument("-Y", help="Y size", type=float, default=None);
    parser.add_argument("-Z", help="Z size", type=float, default=None);
    parser.add_argument("--num-samples",
            help="number of samples along each dimention", type=int, default=2);
    parser.add_argument("--num-samples-X", help="X samples", type=int, default=None);
    parser.add_argument("--num-samples-Y", help="Y samples", type=int, default=None);
    parser.add_argument("--num-samples-Z", help="Z samples", type=int, default=None);
    parser.add_argument("--subdiv", help="subdivision order", type=int,
            default=0);
    parser.add_argument("output", help="output_mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    side_lengths = np.ones(args.dim) * args.size;
    if (args.X is not None):
        side_lengths[0] = args.X
    if (args.Y is not None):
        side_lengths[1] = args.Y
    if (args.dim == 3 and args.Z is not None):
        side_lengths[2] = args.Z

    num_samples = np.ones(args.dim, dtype=int) * args.num_samples;
    if args.num_samples_X is not None:
        num_samples[0] = args.num_samples_X;
    if args.num_samples_Y is not None:
        num_samples[1] = args.num_samples_Y;
    if args.dim == 3 and args.num_samples_Z is not None:
        num_samples[2] = args.num_samples_Z;
    generate_box(side_lengths, num_samples, args.output,
            args.symmetric, args.subdiv);

if __name__ == "__main__":
    main();

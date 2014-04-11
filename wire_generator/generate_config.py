#!/usr/bin/env python

import argparse
import json
from math import sqrt
import numpy as np
from numpy.linalg import solve
import os.path

from WireNetwork import WireNetwork
from WirePattern import WirePattern

def load_wire_network(wire_file):
    network = WireNetwork();
    network.load_from_file(wire_file);
    return network;

def compute_different_types_of_wire_length(wires):
    pattern = WirePattern();
    pattern.set_single_cell_from_wire_network(wires);

    pattern.tile([1, 1, 1]);
    total_length_x1 = pattern.wire_network.total_wire_length;
    pattern.tile([2, 2, 2]);
    total_length_x2 = pattern.wire_network.total_wire_length;
    pattern.tile([3, 3, 3]);
    total_length_x3 = pattern.wire_network.total_wire_length;

    A = np.array(
            [[ 1,  1,   1],
             [ 8,  6, 4.5],
             [27, 18,  12]], dtype=float);
    sol = solve(A, [total_length_x1, total_length_x2, total_length_x3]);

    interior_wire_len = sol[0];
    face_wire_len = sol[1];
    edge_wire_len = sol[2];
    return interior_wire_len, face_wire_len, edge_wire_len;

def compute_cell_size(wires, volume_fraction, thickness):
    interior_wire_len, face_wire_len, edge_wire_len = \
            compute_different_types_of_wire_length(wires);
    wire_length_per_cell = interior_wire_len +\
            0.5 * face_wire_len + 0.25 * edge_wire_len;
    bbox_min, bbox_max = wires.bbox;
    volume = np.prod(bbox_max - bbox_min);
    cell_size = sqrt(wire_length_per_cell * thickness**2 /
            (volume_fraction * volume));

    return (bbox_max - bbox_min) * cell_size;

def compute_repetitions(cell_size, cube_size):
    return np.around(np.divide(cube_size, cell_size)).astype(int).tolist();

def generate_config(wire_file, repeats, thickness, cube_size, output_file):
    """ syntax:
    {
        "wire_network": single_cell_wire_network,
        "thickness": float,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "repeats": [x_reps, y_reps, z_reps],
        "output": output_file
    }
    """
    basename, ext = os.path.splitext(output_file);
    config = {
            "wire_network": wire_file,
            "thickness": thickness,
            "bbox_min": [0.0, 0.0, 0.0],
            "bbox_max": [cube_size, cube_size, cube_size],
            "repeats": repeats,
            "output": basename + ".obj" };
    with open(output_file, 'w') as fout:
        json.dump(config, fout, indent=4)

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate config file for tiling a cube");
    parser.add_argument("--thickness", "-t", help="thickness of each wire",
            default=1.0, type=float);
    parser.add_argument("--size", "-s", help="size of the cube", default=30.0,
            type=float);
    parser.add_argument("--volume-fraction", help="desired volume fraction",
            default=0.25, type=float);
    parser.add_argument("wire_file", help="wire network file");
    parser.add_argument("output_file", help="output config file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();

    wires = load_wire_network(args.wire_file);
    cell_size = compute_cell_size(wires, args.volume_fraction, args.thickness);
    repeats = compute_repetitions(cell_size, args.size);
    generate_config(args.wire_file, repeats, args.thickness, args.size,
            args.output_file);

if __name__ == "__main__":
    main();


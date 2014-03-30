#!/usr/bin/env python
import argparse
import json
import numpy as np
import os.path

from WireInflator import WireInflator
from WireNetwork import WireNetwork
from WirePattern import WirePattern
from timethis import timethis

def parse_config_file(config_file):
    config_dir = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);
        output = config["output"];
        if not os.path.isabs(output):
            config["output"] = os.path.join(config_dir, output);
        return config;

def generate_pattern():
    vertices = [
            # Bottom
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0],
            [-1.0,  0.0,  0.0],
            [ 0.0,  0.0,  1.0],
            [ 0.0,  0.0, -1.0],
            # Top
            [ 0.0,  2.0,  0.0],
            [ 1.0,  2.0,  0.0],
            [-1.0,  2.0,  0.0],
            [ 0.0,  2.0,  1.0],
            [ 0.0,  2.0, -1.0],
            # Middle
            [ 0.5,  1.0,  0.0],
            [ 1.5,  1.0,  0.0],
            [-0.5,  1.0,  0.0],
            [-1.5,  1.0,  0.0],
            [ 0.0,  1.0,  0.5],
            [ 0.0,  1.0,  1.5],
            [ 0.0,  1.0, -0.5],
            [ 0.0,  1.0, -1.5]
            ];
    edges = [
            # Bottom
            [0, 1], [0, 2], [0, 3], [0, 4],
            # Top
            [5, 6], [5, 7], [5, 8], [5, 9],
            # Middle
            [10, 11], [12, 13], [14, 15], [16, 17],
            # Bottom to middle
            [1, 10], [2, 12], [3, 14], [4, 16],
            # Top to middle
            [6, 10], [7, 12], [8, 14], [9, 16]
            ];

    return WirePattern(vertices, edges);

@timethis
def generate_negative_poisson(config):
    """ syntax:
    {
        "thickness": float,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "repeats": [x_reps, y_reps, z_reps],
        "output": output_mesh
    }
    """
    thickness = config.get("thickness");
    bbox_min = np.array(config.get("bbox_min", [0.0, 0.0, 0.0]));
    bbox_max = np.array(config.get("bbox_max", [1.0, 1.0, 1.0]));
    reps = config.get("repeats", [10, 10, 10]);

    pattern = generate_pattern();
    pattern.tile(reps);

    network = pattern.wire_network;
    network_bbox_min, network_bbox_max = network.bbox;
    factor = np.divide(bbox_max - bbox_min, network_bbox_max - network_bbox_min);
    network.scale(factor);

    inflator = WireInflator(network);
    inflator.inflate(thickness);
    inflator.save(str(config["output"]));

def parse_args():
    parser = argparse.ArgumentParser(description="Generate negative Poisson meshes");
    parser.add_argument("config_file", help="configuration file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    config = parse_config_file(args.config_file);
    generate_negative_poisson(config);

if __name__ == "__main__":
    main();
    timethis.summarize();


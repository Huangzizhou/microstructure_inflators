#!/usr/bin/env python
import argparse
import json
import numpy as np
import os.path

from WireInflator import WireInflator
from WireNetwork import WireNetwork

def parse_config_file(config_file):
    config_dir = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);
        output = config["output"];
        if not os.path.isabs(output):
            config["output"] = os.path.join(config_dir, output);
        return config;

def generate_lattice(config):
    """ syntax:
    {
        "thickness": float,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "num_samples": [num_x_samples, num_y_samples, num_z_samples],
        "output": output_mesh
    }
    """
    thickness = config.get("thickness");
    bbox_min = np.array(config.get("bbox_min", [0.0, 0.0, 0.0]));
    bbox_max = np.array(config.get("bbox_max", [1.0, 1.0, 1.0]));
    num_samples = np.array(config.get("num_samples", [10, 10, 10]));
    step_size = np.divide(bbox_max - bbox_min, num_samples-1);

    eps = 1e-3;
    vertices = [];
    for x in np.arange(bbox_min[0], bbox_max[0]+eps, step_size[0]):
        for y in np.arange(bbox_min[1], bbox_max[1]+eps, step_size[1]):
            for z in np.arange(bbox_min[2], bbox_max[2]+eps, step_size[2]):
                vertices.append([x, y, z]);

    edges = [];
    index = lambda i,j,k: i*num_samples[1]*num_samples[2]+j*num_samples[2]+k;
    for i in range(num_samples[0]):
        for j in range(num_samples[1]):
            for k in range(num_samples[2]):
                curr_idx = index(i,j,k);
                if i+1 < num_samples[0]:
                    next_idx_x = index(i+1,j,k);
                    edges.append([curr_idx, next_idx_x]);
                if j+1 < num_samples[1]:
                    next_idx_y = index(i,j+1,k);
                    edges.append([curr_idx, next_idx_y]);
                if k+1 < num_samples[2]:
                    next_idx_z = index(i,j,k+1);
                    edges.append([curr_idx, next_idx_z]);

    wire_network = WireNetwork();
    wire_network.load(vertices, edges);
    inflator = WireInflator(wire_network);
    inflator.inflate(thickness);
    inflator.save(str(config["output"]));

def parse_args():
    parser = argparse.ArgumentParser(description="Generate lattice mesh");
    parser.add_argument("config_file", help="configuration file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    config = parse_config_file(args.config_file);
    generate_lattice(config);

if __name__ == "__main__":
    main();

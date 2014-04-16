#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

from WireNetwork import WireNetwork
from WirePattern import WirePattern
from WireInflator import WireInflator

def parse_config_file(config_file):
    """ syntax:
    {
        "wire_network": single_cell_wire_network,
        "thickness": float,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "repeats": [x_reps, y_reps, z_reps],
        "x_tile_dir": [x, y, z],
        "y_tile_dir": [x, y, z],
        "z_tile_dir": [x, y, z],
        "trim": bool,
        "output": output_file
    }
    """
    config_dir = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);

    def convert_to_abs_path(field_name):
        field = config[field_name];
        if not os.path.isabs(field):
            field = os.path.join(config_dir, field);
        config[field_name] = field;
        
    convert_to_abs_path("wire_network");
    convert_to_abs_path("output");
    return config;

def tile(config):
    network = WireNetwork();
    network.load_from_file(config["wire_network"]);

    pattern = WirePattern();
    pattern.set_single_cell_from_wire_network(network);
    pattern.x_tile_dir = np.array(config.get("x_tile_dir", [1.0, 0.0, 0.0]));
    pattern.y_tile_dir = np.array(config.get("y_tile_dir", [0.0, 1.0, 0.0]));
    pattern.z_tile_dir = np.array(config.get("z_tile_dir", [0.0, 0.0, 1.0]));
    pattern.tile(config["repeats"]);

    tiled_network = pattern.wire_network;
    if config.get("trim", False):
        tiled_network.trim();

    bbox_min, bbox_max = tiled_network.bbox;
    target_bbox_min = np.array(config["bbox_min"]);
    target_bbox_max = np.array(config["bbox_max"]);
    factor = np.divide(target_bbox_max - target_bbox_min, bbox_max - bbox_min);
    tiled_network.scale(factor);

    inflator = WireInflator(tiled_network);
    inflator.inflate(config["thickness"]);
    inflator.save(str(config["output"]));

def parse_args():
    parser = argparse.ArgumentParser(description="Tile a given pattern");
    parser.add_argument("--output", "-o", required=False);
    parser.add_argument("config_file", help="pattern configuration file.");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    config = parse_config_file(args.config_file);
    if args.output is not None:
        config["output"] = args.output;
    tile(config);

if __name__ == "__main__":
    main();


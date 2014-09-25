#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

from core.WireNetwork import WireNetwork
from inflator.WireTiler import WireTiler
from inflator.WireInflator import WireInflator
from inflator.PeriodicWireInflator import PeriodicWireInflator
from parameter.ParameterFactory import ParameterFactory
from utils.find_file import find_file
from utils.timethis import timethis
import PyMesh

def parse_config_file(config_file):
    """ syntax:
    {
        "wire_network": single_cell_wire_network,
        "thickness": float,
        "modifier_file": modifier_file,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "repeats": [x_reps, y_reps, z_reps],
        "hex_mesh": hex_mesh,
        "no_tile": bool,
        "trim": bool,
        "periodic": bool,
        "output": output_file
    }
    """
    config_dir = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);

    def convert_to_abs_path(field_name):
        field = config[field_name];
        if isinstance(field, (unicode, str)):
            config[field_name] = find_file(field, config_dir);

    convert_to_abs_path("wire_network");
    if "hex_mesh" in config:
        convert_to_abs_path("hex_mesh");
    if "modifier_file" in config:
        convert_to_abs_path("modifier_file");
    return config;

def load_parameters(wire_network, default_thickness, modifier_file):
    factory = ParameterFactory(wire_network, default_thickness);
    factory.create_parameters_from_file(modifier_file);
    return factory.parameters;

def tile(config):
    network = load_wire(str(config["wire_network"]));
    parameters = load_parameters(network,
            config["thickness"], config.get("modifier_file"));

    if config.get("no_tile", False):
        tiled_network = network;
    elif "hex_mesh" in config:
        tiled_network = tile_hex(config, network, parameters);
    else:
        tiled_network = tile_box(config, network, parameters);

    if config.get("trim", False):
        tiled_network.trim();

    inflate_and_save(tiled_network, config.get("periodic", False),
            str(config["output"]));

def load_mesh(mesh_file):
    factory = PyMesh.MeshFactory();
    factory.load_file(mesh_file);
    mesh = factory.create();
    return mesh;

def load_wire(wire_file):
    network = WireNetwork();
    network.load_from_file(wire_file);
    network.compute_symmetry_orbits();
    return network;

def inflate_and_save(tiled_network, periodic, output_file):
    if not periodic:
        inflator = WireInflator(tiled_network);
    else:
        inflator = PeriodicWireInflator(tiled_network);

    inflator.inflate();
    inflator.save(output_file);

def tile_hex(config, network, parameters):
    tiler = WireTiler();
    tiler.set_single_cell_from_wire_network(network);
    hex_mesh = load_mesh(str(config["hex_mesh"]));
    tiler.tile_hex_mesh(hex_mesh, parameters);
    tiled_network = tiler.wire_network;
    return tiled_network;

def tile_box(config, network, parameters):
    dim = network.dim;
    tiler = WireTiler();
    tiler.set_single_cell_from_wire_network(network);
    tiler.x_tile_dir = np.array(config.get("x_tile_dir", [1.0, 0.0, 0.0]))[:dim];
    tiler.y_tile_dir = np.array(config.get("y_tile_dir", [0.0, 1.0, 0.0]))[:dim];
    tiler.z_tile_dir = np.array(config.get("z_tile_dir", [0.0, 0.0, 1.0]))[:dim];
    tiler.tile(config["repeats"][:dim], parameters);

    tiled_network = tiler.wire_network;

    bbox_min, bbox_max = tiled_network.bbox;
    bbox_size = bbox_max - bbox_min;
    non_zero_dim = bbox_size > 0.0;

    target_bbox_min = np.array(config["bbox_min"])[:dim];
    target_bbox_max = np.array(config["bbox_max"])[:dim];
    target_bbox_size = target_bbox_max - target_bbox_min;

    factor = np.ones(dim);
    factor[non_zero_dim] = np.divide(
            target_bbox_size[non_zero_dim], bbox_size[non_zero_dim]);
    tiled_network.scale(factor);
    return tiled_network;

def parse_args():
    parser = argparse.ArgumentParser(description="Tile a given pattern");
    parser.add_argument("--output", "-o", required=True);
    parser.add_argument("--timing", help="display running times",
            action="store_true");
    parser.add_argument("config_file", help="pattern configuration file.");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    config = parse_config_file(args.config_file);
    if args.output is not None:
        config["output"] = args.output;
    tile(config);
    if args.timing:
        timethis.summarize();

if __name__ == "__main__":
    main();


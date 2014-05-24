#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

import PyMeshSetting
import PyMesh

from WireNetwork import WireNetwork
from WirePattern import WirePattern
from WireInflator import WireInflator
from WireModifierFactory import WireModifierFactory
from PeriodicWireInflator import PeriodicWireInflator
from timethis import timethis

def parse_config_file(config_file):
    """ syntax:
    {
        "wire_network": single_cell_wire_network,
        "thickness": float,
        "modifiers": modifier_file,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "repeats": [x_reps, y_reps, z_reps],
        "x_tile_dir": [x, y, z],
        "y_tile_dir": [x, y, z],
        "z_tile_dir": [x, y, z],
        "hex_mesh": hex_mesh,
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
            if not os.path.isabs(field):
                field = os.path.join(config_dir, field);
            config[field_name] = field;

    convert_to_abs_path("wire_network");
    convert_to_abs_path("output");
    if "vertex_offset" in config:
        convert_to_abs_path("vertex_offset");
    if "hex_mesh" in config:
        convert_to_abs_path("hex_mesh");
    if "modifier_file" in config:
        convert_to_abs_path("modifier_file");
    return config;

def set_uniform_thickness(wire_network, thickness):
    thickness = np.ones(wire_network.num_vertices) * thickness;
    wire_network.attributes.add("vertex_thickness", thickness);

def load_modifiers(modifier_file):
    if modifier_file is None:
        return [];

    modifiers = WireModifierFactory.create_from_file(str(modifier_file));
    return modifiers;

def tile(config):
    network = load_wire(str(config["wire_network"]));
    set_uniform_thickness(network, config["thickness"]);
    modifiers = load_modifiers(config.get("modifier_file", None));

    if "hex_mesh" in config:
        tiled_network = tile_hex(config, network, modifiers);
    else:
        tiled_network = tile_box(config, network, modifiers);

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
    return network;

def inflate_and_save(tiled_network, periodic, output_file):
    if not periodic:
        inflator = WireInflator(tiled_network);
    else:
        inflator = PeriodicWireInflator(tiled_network);

    inflator.inflate();
    inflator.save(output_file);

def tile_hex(config, network, modifiers):
    pattern = WirePattern();
    pattern.set_single_cell_from_wire_network(network);
    hex_mesh = load_mesh(str(config["hex_mesh"]));
    pattern.tile_hex_mesh(hex_mesh, modifiers);

    tiled_network = pattern.wire_network;
    return tiled_network;

def tile_box(config, network, modifiers):
    pattern = WirePattern();
    pattern.set_single_cell_from_wire_network(network);
    pattern.x_tile_dir = np.array(config.get("x_tile_dir", [1.0, 0.0, 0.0]));
    pattern.y_tile_dir = np.array(config.get("y_tile_dir", [0.0, 1.0, 0.0]));
    pattern.z_tile_dir = np.array(config.get("z_tile_dir", [0.0, 0.0, 1.0]));
    pattern.tile(config["repeats"], modifiers);

    tiled_network = pattern.wire_network;

    bbox_min, bbox_max = tiled_network.bbox;
    target_bbox_min = np.array(config["bbox_min"]);
    target_bbox_max = np.array(config["bbox_max"]);
    factor = np.divide(target_bbox_max - target_bbox_min, bbox_max - bbox_min);
    tiled_network.scale(factor);
    return tiled_network;

def parse_args():
    parser = argparse.ArgumentParser(description="Tile a given pattern");
    parser.add_argument("--output", "-o", required=False);
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


#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path
import os
import sys

sys.path.append("../lib");
sys.path.append("../swig");

import PyWireInflator2D
import PyMeshSetting
import PyMesh

from WireModifierFactory import WireModifierFactory

def form_mesh(vertices, faces, voxels=np.array([])):
    dim = vertices.shape[1];
    factory = PyMesh.MeshFactory();
    if dim == 3:
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels.ravel(order="C"), 3, 3, 4);
    elif dim == 2:
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels.ravel(order="C"), 2, 3, 4);
    return factory.create();

def save_mesh(filename, mesh, *attributes):
    writer = PyMesh.MeshWriter.create_writer(filename);
    for attr in attributes:
        if not mesh.has_attribute(attr):
            raise KeyError("Attribute {} is not found in mesh".format(attr));
        writer.with_attribute(attr);
    writer.write_mesh(mesh);

def load_mesh(filename):
    factory = PyMesh.MeshFactory();
    factory.load_file(filename);
    return factory.create();

def parse_config_file(config_file):
    """ syntax:
    {
        "wire_network": single_cell_wire_network,
        "thickness": float,
        "modifier_file": modifier_file,
        "bbox_min": [min_x, min_y],
        "bbox_max": [max_x, max_y],
        "repeats": [x_reps, y_reps],
        "quad_mesh": quad_mesh,
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
    modifiers = load_modifiers(config.get("modifier_file", None));
    if "quad_mesh" in config:
        rows, cols, param = tile_quad(config, modifiers);
    else:
        rows, cols, param, scale_factor = tile_box(config, modifiers);

    if config.get("trim", False):
        raise NotImplementedError("Trimming is not supported");

    inflator = PyWireInflator2D.WireInflatorFacade(rows, cols);
    for i in range(rows):
        for j in range(cols):
            p = param[i][j];
            if p is not None:
                inflator.set_parameter(i,j,p);

    if config.get("periodic", False):
        inflator.generate_periodic_pattern();
    else:
        inflator.generate_tiled_pattern();

    vertices = inflator.get_vertices().reshape((-1, 2), order="C");
    vertices = vertices * scale_factor;
    faces = inflator.get_triangles().reshape((-1, 3), order="C");

    mesh = form_mesh(vertices, faces, np.array([]));
    save_mesh(str(config["output"]), mesh);

def tile_quad(config, modifiers):
    raise NotImplementedError("Tiling quad is not implemented yet");
    #pattern = WirePattern();
    #pattern.set_single_cell_from_wire_network(network);
    #hex_mesh = load_mesh(str(config["hex_mesh"]));
    #pattern.tile_hex_mesh(hex_mesh, modifiers);

    #tiled_network = pattern.wire_network;
    #return tiled_network;

def tile_box(config, modifiers):
    rows, cols = config["repeats"];
    bbox_min = np.array(config["bbox_min"]);
    bbox_max = np.array(config["bbox_max"]);
    scale_factor = np.divide(bbox_max - bbox_min, [rows, cols]);

    default_parameter = np.zeros(8);
    default_parameter[:5] = config["thickness"];

    for modifier in modifiers:
        modifier.modify(default_parameter);
        print(default_parameter);

    params = [[default_parameter for i in range(cols)] for j in range(rows)];
    return rows, cols, params, scale_factor;

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


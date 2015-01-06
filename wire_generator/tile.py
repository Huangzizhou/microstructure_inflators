#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

from core.WireNetwork import WireNetwork
from inflator.InflatorFacade import InflatorFacade
from parameter.PyParameters import PyParameters
from utils.find_file import find_file
from utils.timethis import timethis
import PyMesh

def parse_config_file(config_file):
    """ syntax:
    {
        # opitional options:
        "trim": bool,
        "periodic": bool,
        "subdiv": #,
        "subdiv_method": "simple" or "loop"
        "rel_geometry_correction": [#, #, #],
        "abs_geometry_correction": [#, #, #],
        "geometry_correction_cap": #,

        # options for specifying parameters
        "thickness": float,
        "modifier_file": modifier_file,
        "dof_file": dof_file,

        # options needed by guide bbox
        "wire_network": single_cell_wire_network,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "repeats": [x_reps, y_reps, z_reps],

        # options needed by guide mesh
        "wire_network": single_cell_wire_network,
        "guide_mesh": guide_mesh,
        "dof_type": "isotropic" | "orthotropic",
        "thickness_type": "vertex" | "edge"

        # options neede by mixed pattern tiling
        "guide_mesh": guide_mesh,
        "wire_list_file": wire_list_file,
    }
    """
    config_dir = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);

    def convert_to_abs_path(field_name):
        field = config[field_name];
        if isinstance(field, (unicode, str)):
            config[field_name] = str(find_file(field, config_dir));

    if "wire_network" in config:
        convert_to_abs_path("wire_network");
    if "guide_mesh" in config:
        convert_to_abs_path("guide_mesh");
    if "modifier_file" in config:
        convert_to_abs_path("modifier_file");
    if "dof_file" in config:
        convert_to_abs_path("dof_file");
    if "wire_list_file" in config:
        convert_to_abs_path("wire_list_file");
    return config;

def load_mesh(mesh_file):
    factory = PyMesh.MeshFactory();
    factory.load_file(str(mesh_file));
    factory.drop_zero_dim();
    mesh = factory.create();
    return mesh;

def save_mesh(mesh, mesh_file):
    basename, ext = os.path.splitext(mesh_file);
    writer = PyMesh.MeshWriter.create_writer(mesh_file);
    if ext in (".msh", ".ply"):
        attribute_names = mesh.get_attribute_names();
        for attr_name in attribute_names:
            writer.with_attribute(attr_name);

    writer.write_mesh(mesh);

def load_wire(wire_file):
    network = WireNetwork();
    network.load_from_file(wire_file);
    network.compute_symmetry_orbits();
    return network;

def load_wires(wire_list_file):
    root_dir = os.path.dirname(wire_list_file);
    with open(wire_list_file, 'r') as fin:
        wire_files = [str(name.strip()) for name in fin];
    wire_files = [find_file(name, root_dir) for name in wire_files];
    wires = [load_wire(name) for name in wire_files];
    return wires;

def load_parameters(wire_network, config):
    parameters = PyParameters(wire_network, config["thickness"]);
    if "modifier_file" in config:
        if "dof_file" in config:
            print("dof_file is shadowed by modifier file!");
        parameters.load_modifier_file(config["modifier_file"]);
    elif "dof_file" in config:
        parameters.load_dof_file(config["dof_file"]);
    return parameters;

def extract_options(config):
    options = {
            "trim": config.get("trim", False),
            "periodic": config.get("periodic", False),
            "subdiv": config.get("subdiv", 1),
            "subdiv_method": str(config.get("subdiv_method", "simple")),
            "rel_geometry_correction": config.get("rel_geometry_correction"),
            "abs_geometry_correction": config.get("abs_geometry_correction"),
            "geometry_correction_cap": config.get("geometry_correction_cap"),
            }
    return options;

@timethis
def tile(config):
    if "guide_mesh" in config:
        if "wire_list_file" in config:
            return tile_with_mixed_patterns(config);
        else:
            return tile_with_guide_mesh(config);
    else:
        return tile_with_guide_box(config);

def tile_with_guide_box(config):
    options = extract_options(config);
    network = load_wire(str(config["wire_network"]));
    parameters = load_parameters(network, config);
    inflator_driver = InflatorFacade.create(network, parameters);
    mesh = inflator_driver.inflate_with_guide_box(
            config["bbox_min"], config["bbox_max"],
            config["repeats"], options);
    return mesh;

def tile_with_guide_mesh(config):
    options = extract_options(config);
    network = load_wire(str(config["wire_network"]));
    parameters = load_parameters(network, config);
    inflator_driver = InflatorFacade.create(network, parameters);
    guide_mesh = load_mesh(config["guide_mesh"]);
    mesh = inflator_driver.inflate_with_guide_mesh(guide_mesh, options);
    return mesh;

def tile_with_mixed_patterns(config):
    options = extract_options(config);
    options["dof_type"] = str(config.get("dof_type", "isotropic"));
    options["thickness_type"] = str(config.get("thickness_type", "vertex"));
    networks = load_wires(str(config["wire_list_file"]));
    guide_mesh = load_mesh(config["guide_mesh"]);
    inflator_driver = InflatorFacade.create_mixed(networks);
    mesh = inflator_driver.inflate_with_mixed_patterns(guide_mesh, options);
    return mesh;

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
    mesh = tile(config);
    save_mesh(mesh, config["output"]);
    if args.timing:
        timethis.summarize();

if __name__ == "__main__":
    main();


#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

from core.WireNetwork import WireNetwork
from inflator.InflatorFacade import InflatorFacade
#from parameter.ParameterFactory import ParameterFactory
from parameter.PyParameters import PyParameters
from utils.find_file import find_file
from utils.timethis import timethis
import PyMesh

def parse_config_file(config_file):
    """ syntax:
    {
        "wire_network": single_cell_wire_network,
        "thickness": float,
        "modifier_file": modifier_file,
        "dof_file": dof_file,
        "bbox_min": [min_x, min_y, min_z],
        "bbox_max": [max_x, max_y, max_z],
        "repeats": [x_reps, y_reps, z_reps],
        "guide_mesh": guide_mesh,
        "no_tile": bool,
        "trim": bool,
        "periodic": bool,
        "subdiv": #,
        "subdiv_method": "simple" or "loop"
        "geometry_correction": [#, #, #],
        "output": output_file
    }
    """
    config_dir = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);

    def convert_to_abs_path(field_name):
        field = config[field_name];
        if isinstance(field, (unicode, str)):
            config[field_name] = str(find_file(field, config_dir));

    convert_to_abs_path("wire_network");
    if "guide_mesh" in config:
        convert_to_abs_path("guide_mesh");
    if "modifier_file" in config:
        convert_to_abs_path("modifier_file");
    if "dof_file" in config:
        convert_to_abs_path("dof_file");
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

def load_parameters_old(wire_network, default_thickness, modifier_file):
    factory = ParameterFactory(wire_network, default_thickness);
    factory.create_parameters_from_file(modifier_file);
    return factory.parameters;

def load_parameters(wire_network, config):
    parameters = PyParameters(wire_network, config["thickness"]);
    if "modifier_file" in config:
        if "dof_file" in config:
            print("dof_file is shadowed by modifier file!");
        parameters.load_modifier_file(config["modifier_file"]);
    elif "dof_file" in config:
        parameters.load_dof_file(config["dof_file"]);
    return parameters;

@timethis
def tile(config):
    network = load_wire(str(config["wire_network"]));
    parameters = load_parameters(network, config);

    options = {
            "trim": config.get("trim", False),
            "periodic": config.get("periodic", False),
            "subdiv": config.get("subdiv", 1),
            "subdiv_method": str(config.get("subdiv_method", "simple")),
            "geometry_correction": config.get(
                "geometry_correction", np.zeros(network.dim)),
            }
    inflator_driver = InflatorFacade.create(network, parameters);
    if "guide_mesh" in config:
        guide_mesh = load_mesh(config["guide_mesh"]);
        mesh = inflator_driver.inflate_with_guide_mesh(guide_mesh, options);
    else:
        mesh = inflator_driver.inflate_with_guide_box(
                config["bbox_min"], config["bbox_max"],
                config["repeats"], options);
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


#!/usr/bin/env python

import argparse
import json
import numpy as np

from core.WireNetwork import WireNetwork
import PyMesh
import PyWires

def load_wire(wire_file):
    network = WireNetwork();
    network.load_from_file(wire_file);
    network.compute_symmetry_orbits();
    return network;

def generate_vertex_thickness_config(wire_network, thickness_dof_map,
        default_thickness):
    vertex_orbits = wire_network.attributes["orthotropic_symmetry_vertex_orbit"];
    vertex_orbits = vertex_orbits.astype(int);
    orbit_indices = np.unique(vertex_orbits).tolist();
    orbit_dof = [];
    for idx in orbit_indices:
        dof_indices = thickness_dof_map[vertex_orbits == idx];
        dof_id = dof_indices[0];
        if not np.all(dof_indices == dof_id):
            raise RuntimeError(
                    "Vertex symmetry orbit incompatible with C++ version");
        if dof_id == -1:
            orbit_dof.append(default_thickness);
        else:
            orbit_dof.append("{{dof_{}}}".format(dof_id));

    config = {
            "type": "vertex_orbit",
            "default": default_thickness,
            "effective_orbits": orbit_indices,
            "thickness": orbit_dof
            };
    return config;

def generate_edge_thickness_config(wire_network, thickness_dof_map,
        default_thickness):
    edge_orbits = wire_network.attributes["orthotropic_symmetry_edge_orbit"];
    edge_orbits = edge_orbits.astype(int);
    orbit_indices = np.unique(edge_orbits).tolist();
    orbit_dof = [];
    for idx in orbit_indices:
        dof_indices = thickness_dof_map[edge_orbits == idx];
        dof_id = dof_indices[0];
        if not np.all(dof_indices == dof_id):
            raise RuntimeError(
                    "Edge symmetry orbit incompatible with C++ version");
        if dof_id == -1:
            orbit_dof.append(default_thickness);
        else:
            orbit_dof.append("{{dof_{}}}".format(dof_id));

    config = {
            "type": "edge_orbit",
            "default": default_thickness,
            "effective_orbits": orbit_indices,
            "thickness": orbit_dof
            };
    return config;

def generate_vertex_offset_config(wire_network, offset_dof_map):
    dim = wire_network.dim;
    vertex_orbits = wire_network.attributes["orthotropic_symmetry_vertex_orbit"];
    vertex_orbits = vertex_orbits.astype(int);
    orbit_indices = np.unique(vertex_orbits).tolist();

    orbit_dof = [];
    for idx in orbit_indices:
        orbit_offset = [];
        for axis in range(dim):
            dof_indices = offset_dof_map[vertex_orbits == idx, axis].ravel();
            dof_id = dof_indices[0];
            if not np.all(dof_indices == dof_id):
                raise RuntimeError(
                        "Vertex symmetry orbit incompatible with C++ version");

            if dof_id == -1:
                orbit_offset.append(0.0);
            else:
                orbit_offset.append("{{dof_{}}}".format(dof_id));
        orbit_dof.append(orbit_offset);

    config = {
            "type": "vertex_orbit",
            "effector_orbits": orbit_indices,
            "offset_percentages": orbit_dof
            };
    return config;

def extract_parameters(wire_network, default_thickness):
    c_parameters = PyWires.ParameterManager.create(
            wire_network.raw_wires, default_thickness);
    thickness_dof_map = c_parameters.get_thickness_dof_map().ravel();
    offset_dof_map = c_parameters.get_offset_dof_map();

    config = { "orbit_type": "orthotropic" };
    if c_parameters.get_thickness_type() == PyWires.VERTEX:
        thickness_config = generate_vertex_thickness_config(
                        wire_network, thickness_dof_map, default_thickness);
    else:
        thickness_config = generate_edge_thickness_config(
                wire_network, thickness_dof_map, default_thickness);
    config["thickness"] = thickness_config;
    config["vertex_offset"] = generate_vertex_offset_config(
            wire_network, offset_dof_map);
    return config;

def dump_json(filename, config):
    with open(filename, 'w') as fout:
        json.dump(config, fout, indent=4);

def parse_args():
    parser = argparse.ArgumentParser(
            description="Automatically extract all parameters from a wire mesh");
    parser.add_argument("--default-thickness", default=0.5, type=float);
    parser.add_argument("wire_file", help="input wire file");
    parser.add_argument("modifier_file", help="output modifier_file");
    return parser.parse_args();

def main():
    args = parse_args();

    wire_network = load_wire(args.wire_file);
    config = extract_parameters(wire_network, args.default_thickness);
    dump_json(args.modifier_file, config);

if __name__ == "__main__":
    main();


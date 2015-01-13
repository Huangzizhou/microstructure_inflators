#!/usr/bin/env python

import argparse
import json
import os.path

from core.WireNetwork import WireNetwork

def load(wire_file):
    wire_network = WireNetwork();
    wire_network.load_from_file(wire_file);
    return wire_network;

def extract_orbits(wire_network):
    isotropic_vertex_attr_name = "vertex_cubic_symmetry_orbit";
    orthotropic_vertex_attr_name = "vertex_symmetry_orbit";
    isotropic_edge_attr_name = "edge_cubic_symmetry_orbit";
    orthotropic_edge_attr_name = "edge_symmetry_orbit";

    wire_network.add_attribute(orthotropic_vertex_attr_name);
    orthotropic_vertex_orbits = wire_network.get_attribute(
            orthotropic_vertex_attr_name).ravel().astype(int).tolist();

    wire_network.add_attribute(isotropic_vertex_attr_name);
    isotropic_vertex_orbits = wire_network.get_attribute(
            isotropic_vertex_attr_name).ravel().astype(int).tolist();

    wire_network.add_attribute(orthotropic_edge_attr_name);
    orthotropic_edge_orbits = wire_network.get_attribute(
            orthotropic_edge_attr_name).ravel().astype(int).tolist();

    wire_network.add_attribute(isotropic_edge_attr_name);
    isotropic_edge_orbits = wire_network.get_attribute(
            isotropic_edge_attr_name).ravel().astype(int).tolist();

    return orthotropic_vertex_orbits, isotropic_vertex_orbits,\
            orthotropic_edge_orbits, isotropic_edge_orbits;

def generate_index_map(orbits):
    orbit_map = {};
    for i,orbit_id in enumerate(orbits):
        if orbit_id in orbit_map:
            orbit_map[orbit_id].append(i);
        else:
            orbit_map[orbit_id] = [i];
    return orbit_map;

def save_orbits(orbit_file,
        orthotropic_vertex_orbits,
        isotropic_vertex_orbits,
        orthotropic_edge_orbits,
        isotropic_edge_orbits):
    orthotropic_vertex_orbit_map = generate_index_map(orthotropic_vertex_orbits);
    isotropic_vertex_orbit_map = generate_index_map(isotropic_vertex_orbits);
    orthotropic_edge_orbit_map = generate_index_map(orthotropic_edge_orbits);
    isotropic_edge_orbit_map = generate_index_map(isotropic_edge_orbits);

    contents = {
            "orthotropic_vertex_orbits": orthotropic_vertex_orbit_map.values(),
            "isotropic_vertex_orbits": isotropic_vertex_orbit_map.values(),
            "orthotropic_edge_orbits": orthotropic_edge_orbit_map.values(),
            "isotropic_edge_orbits": isotropic_edge_orbit_map.values()
            };
    with open(orbit_file, 'w') as fout:
        json.dump(contents, fout, indent=4);

def parse_args():
    parser = argparse.ArgumentParser(
            description="Extract symmetry orbits from wire files");
    parser.add_argument("wire_files", help="wire files", nargs="+");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    for wire_file in args.wire_files:
        basename, ext = os.path.splitext(wire_file);
        orbit_file = basename + ".orbit";

        wire_network = load(wire_file);
        orbits = extract_orbits(wire_network);
        save_orbits(orbit_file, *orbits);

if __name__ == "__main__":
    main();


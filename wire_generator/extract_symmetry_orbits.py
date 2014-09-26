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
    wire_network.attributes.add("symmetry_vertex_orbit");
    vertex_orbits = wire_network.attributes["symmetry_vertex_orbit"];
    wire_network.attributes.add("symmetry_edge_orbit");
    edge_orbits = wire_network.attributes["symmetry_edge_orbit"];
    return vertex_orbits, edge_orbits;

def generate_index_map(orbits):
    orbit_map = {};
    for i,orbit_id in enumerate(orbits):
        if orbit_id in orbit_map:
            orbit_map[orbit_id].append(i);
        else:
            orbit_map[orbit_id] = [i];
    return orbit_map;

def save_orbits(orbit_file, vertex_orbits, edge_orbits):
    vertex_orbit_map = generate_index_map(vertex_orbits);
    edge_orbit_map = generate_index_map(edge_orbits);

    contents = {
            "vertex_orbits": vertex_orbit_map.values(),
            "edge_orbits": edge_orbit_map.values()
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
        vertex_orbits, edge_orbits = extract_orbits(wire_network);
        save_orbits(orbit_file, vertex_orbits, edge_orbits);

if __name__ == "__main__":
    main();


#!/usr/bin/env python

import argparse
import json
import os.path

from WireNetwork import WireNetwork

def load(wire_file):
    wire_network = WireNetwork();
    wire_network.load_from_file(wire_file);
    return wire_network;

def extract_orbits(wire_network):
    wire_network.attributes.add("symmetry_orbit");
    orbits = wire_network.attributes["symmetry_orbit"];
    return orbits;

def save_orbits(orbit_file, orbits):
    orbit_map = {};
    for i,orbit_id in enumerate(orbits):
        if orbit_id in orbit_map:
            orbit_map[orbit_id].append(i);
        else:
            orbit_map[orbit_id] = [i];

    contents = {"orbits": orbit_map.values()};
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
        save_orbits(orbit_file, orbits);

if __name__ == "__main__":
    main();


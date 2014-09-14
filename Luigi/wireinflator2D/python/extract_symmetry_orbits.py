#!/usr/bin/env python
import argparse
import json
import numpy as np
import os.path

import PyWireInflator2DSetting
import PyWireInflator2D

def save_orbits(wire_file, orbits):
    basename, ext = os.path.splitext(wire_file);
    orbit_file = basename + ".orbit";
    with open(orbit_file, 'w') as fout:
        json.dump({ "vertex_orbits": orbits }, fout);

def parse_args():
    parser = argparse.ArgumentParser(
            description="extract vertex symmetry orbits");
    parser.add_argument("wire_files", nargs="+");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    for wire_file in args.wire_files:
        inflator = PyWireInflator2D.WireInflatorFacade(wire_file);
        num_params = inflator.get_num_parameters();

        orbits = [];
        for i in range(num_params):
            orbit = inflator.get_affected_vertex_orbit(i).ravel().tolist();
            orbits.append(tuple(orbit));

        orbits = np.unique(orbits);
        orbits = [list(orbit) for orbit in orbits]
        save_orbits(wire_file, orbits);

if __name__ == "__main__":
    main();

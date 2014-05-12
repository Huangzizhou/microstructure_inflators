#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

from WireNetwork import WireNetwork

def load(wire_file):
    wire_network = WireNetwork();
    wire_network.load_from_file(wire_file);
    return wire_network;

def extract_orbits(wire_network):
    wire_network.attributes.add("symmetry_orbit");
    orbits = wire_network.attributes["symmetry_orbit"];

    orbit_map = {};
    for i,orbit_id in enumerate(orbits):
        if orbit_id in orbit_map:
            orbit_map[orbit_id].append(i);
        else:
            orbit_map[orbit_id] = [i];
    return orbit_map.values();

def thickness_sweep(wire_network, orbits, min_thickness, max_thickness, num_samples):
    num_vertices = len(wire_network.vertices);
    num_orbits = len(orbits);
    thickness_space = [np.linspace(min_thickness, max_thickness, num_samples)]\
            * num_orbits;
    grid = np.meshgrid(*thickness_space);
    parameters = np.vstack(map(np.ravel, grid)).T;

    per_vertex_thicknesses = [];
    for param in parameters:
        thickness = np.zeros(num_vertices);
        for orbit, t in zip(orbits, param):
            thickness[orbit] = t;
        per_vertex_thicknesses.append(thickness);
    return per_vertex_thicknesses, parameters;

def save_thickness(per_vertex_thicknesses, parameters, param_file):
    basename, ext = os.path.splitext(param_file);
    for i, thickness in enumerate(per_vertex_thicknesses):
        filename = "{}_{}{}".format(basename, i, ext);
        param = parameters[i];
        contents = {
                "thickness": thickness.tolist(),
                "parameter": param.tolist()
                };
        with open(filename, 'w') as fout:
            json.dump(contents, fout, indent=4);

def parse_args():
    parser = argparse.ArgumentParser(
            description="Generate config file specifying thickness and offset");
    parser.add_argument("--min-thickness", help="minimum wire thickness",
            type=float, default=0.3);
    parser.add_argument("--max-thickness", help="maximum wire thickness",
            type=float, default=1.0);
    parser.add_argument("--num-thickness-samples",
            help="number of thickness to sample",
            type=int, default=3);
    parser.add_argument("wire_file", help="single cell wire file");
    parser.add_argument("param_file", help="output parameter file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();

    wire_network = load(args.wire_file);
    orbits = extract_orbits(wire_network);
    per_vertex_thicknesses, parameters = thickness_sweep(wire_network, orbits,
            args.min_thickness, args.max_thickness, args.num_thickness_samples);
    save_thickness(per_vertex_thicknesses, parameters, args.param_file);

if __name__ == "__main__":
    main();

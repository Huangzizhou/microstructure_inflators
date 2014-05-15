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

def extract_orbits(wire_network, per_edge):
    if per_edge:
        wire_network.attributes.add("symmetry_edge_orbit");
        orbits = wire_network.attributes["symmetry_edge_orbit"];
    else:
        wire_network.attributes.add("symmetry_vertex_orbit");
        orbits = wire_network.attributes["symmetry_vertex_orbit"];

    orbit_map = {};
    for i,orbit_id in enumerate(orbits):
        if orbit_id in orbit_map:
            orbit_map[orbit_id].append(i);
        else:
            orbit_map[orbit_id] = [i];
    return orbit_map.values();

def thickness_sweep(num_entries, orbits, min_thickness, max_thickness, num_samples):
    num_orbits = len(orbits);
    thickness_space = [np.linspace(min_thickness, max_thickness, num_samples)]\
            * num_orbits;
    grid = np.meshgrid(*thickness_space);
    parameters = np.vstack(map(np.ravel, grid)).T;

    thicknesses = [];
    for param in parameters:
        thickness = np.zeros(num_entries);
        for orbit, t in zip(orbits, param):
            thickness[orbit] = t;
        thicknesses.append(thickness);
    return thicknesses, parameters;

def save_thickness(thicknesses, parameters, param_file, per_edge):
    if per_edge:
        thickness_key = "edge_thickness";
    else:
        thickness_key = "vertex_thickness";

    basename, ext = os.path.splitext(param_file);
    for i, thickness in enumerate(thicknesses):
        filename = "{}_{}{}".format(basename, i, ext);
        param = parameters[i];
        contents = {
                thickness_key: thickness.tolist(),
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
    parser.add_argument("--per-edge", help="assign thickness value to edge",
            action="store_true");
    parser.add_argument("wire_file", help="single cell wire file");
    parser.add_argument("param_file", help="output parameter file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();

    wire_network = load(args.wire_file);

    if args.per_edge:
        num_params = wire_network.num_edges;
    else:
        num_params = wire_network.num_vertices;

    orbits = extract_orbits(wire_network, args.per_edge);
    thicknesses, parameters = thickness_sweep(num_params, orbits,
            args.min_thickness, args.max_thickness, args.num_thickness_samples);
    save_thickness(thicknesses, parameters, args.param_file, args.per_edge);

if __name__ == "__main__":
    main();

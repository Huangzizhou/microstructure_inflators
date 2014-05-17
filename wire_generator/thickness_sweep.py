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

def generate_index_map(orbits):
    orbit_map = {};
    for i,orbit_id in enumerate(orbits):
        if orbit_id in orbit_map:
            orbit_map[orbit_id].append(i);
        else:
            orbit_map[orbit_id] = [i];
    return orbit_map;

def extract_vertex_orbits(wire_network):
    wire_network.attributes.add("symmetry_vertex_orbit");
    orbits = wire_network.attributes["symmetry_vertex_orbit"];
    orbit_map = generate_index_map(orbits);
    return orbit_map.values();

def extract_edge_orbits(wire_network):
    wire_network.attributes.add("symmetry_edge_orbit");
    orbits = wire_network.attributes["symmetry_edge_orbit"];
    orbit_map = generate_index_map(orbits);
    return orbit_map.values();

def thickness_sweep(orbits, min_thickness, max_thickness, num_samples):
    num_orbits = len(orbits);
    num_entries = np.sum(map(len, orbits));

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

class VertexOrbit:
    def __init__(self, wire_network, orbit):
        self.tol = 1e-6;
        self.wire_network = wire_network;
        self.orbit = orbit;
        self.__compute_orbit_dof_map();

    def __compute_orbit_dof_map(self):
        bbox_min, bbox_max = self.wire_network.bbox;
        orbit_vertices = self.wire_network.vertices[self.orbit];
        self.orbit_bbox_min = np.amin(orbit_vertices, axis=0);
        self.orbit_bbox_max = np.amax(orbit_vertices, axis=0);
        self.orbit_bbox_size = self.orbit_bbox_max - self.orbit_bbox_min;
        self.orbit_bbox_center = 0.5 *\
                (self.orbit_bbox_min + self.orbit_bbox_max);

        zero_dim = np.absolute(self.orbit_bbox_size) < self.tol;
        on_min_boundary = np.absolute(self.orbit_bbox_min - bbox_min) < self.tol;
        on_max_boundary = np.absolute(self.orbit_bbox_max - bbox_max) < self.tol;

        neg_dof_map = np.logical_or(zero_dim, on_min_boundary);
        neg_dof_map = np.logical_or(neg_dof_map, on_max_boundary);
        self.dof_map = np.logical_not(neg_dof_map);

    def get_vertex_offsets(self, percentages):
        assert(len(percentages) == self.num_dof);
        ori_vertices = self.wire_network.vertices[self.orbit];
        offset = (ori_vertices - self.orbit_bbox_center);
        offset[:,self.dof_map] *= percentages;
        offset[:,np.logical_not(self.dof_map)] = 0;
        return offset;

    @property
    def num_dof(self):
        return np.sum(self.dof_map);

def offset_sweep(wire_network, vertex_orbits, max_offset_percent, num_samples):
    num_orbits = len(vertex_orbits)
    orbits = [VertexOrbit(wire_network, orbit) for orbit in vertex_orbits];
    num_dof = np.sum([orbit.num_dof for orbit in orbits]);
    print("offset dof: {}".format(num_dof));

    percentages = np.linspace(-max_offset_percent, max_offset_percent, num_samples);
    offset_percent_space = [percentages] * num_dof;
    grid = np.meshgrid(*offset_percent_space);
    parameters = np.vstack(map(np.ravel, grid)).T;

    offsets = [];
    for param in parameters:
        offset = np.zeros((wire_network.num_vertices, wire_network.dim));
        dof_count = 0;
        for orbit in orbits:
            orbit_num_dof = orbit.num_dof;
            orbit_dof_offset = orbit.get_vertex_offsets(
                    param[dof_count:dof_count+orbit_num_dof]);
            offset[orbit.orbit] = orbit_dof_offset;
            dof_count += orbit_num_dof;
        offsets.append(offset);

    return offsets, parameters;

def save_parameters(param_file, per_edge, thicknesses, thickness_parameters,
        offsets, offset_parameters):
    if per_edge:
        thickness_key = "edge_thickness";
    else:
        thickness_key = "vertex_thickness";

    basename, ext = os.path.splitext(param_file);

    num_thicknesses = len(thicknesses);
    num_offsets = len(offsets);
    print(num_thicknesses);
    print(num_offsets);

    import itertools
    for i,j in itertools.product(range(num_thicknesses), range(num_offsets)):
        thickness = thicknesses[i];
        thickness_param = thickness_parameters[i];
        offset = offsets[j];
        offset_param = offset_parameters[j];

        contents = {
                thickness_key: thickness.tolist(),
                "thickness_parameter": thickness_param.tolist(),
                "vertex_offset": offset.tolist(),
                "offset_parameter": offset_param.tolist(),
                };

        filename = "{}_thickness_{}_offset_{}{}".format(basename, i, j, ext);
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
    parser.add_argument("--max-offset-percent",
            help="maximum offset percentage",
            type=float, default=0.2);
    parser.add_argument("--num-offset-samples",
            help="number of offsets to samples", type=int, default=3);
    parser.add_argument("--per-edge", help="assign thickness value to edge",
            action="store_true");
    parser.add_argument("wire_file", help="single cell wire file");
    parser.add_argument("param_file", help="output parameter file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();

    wire_network = load(args.wire_file);

    vertex_orbits = extract_vertex_orbits(wire_network);
    edge_orbits = extract_edge_orbits(wire_network);

    orbits = edge_orbits if args.per_edge else vertex_orbits;

    thicknesses, thickness_parameters = thickness_sweep(orbits,
            args.min_thickness, args.max_thickness, args.num_thickness_samples);

    offsets, offset_parameters = offset_sweep(wire_network, vertex_orbits,
            args.max_offset_percent, args.num_offset_samples);

    save_parameters(
            args.param_file, args.per_edge,
            thicknesses, thickness_parameters,
            offsets, offset_parameters);

if __name__ == "__main__":
    main();

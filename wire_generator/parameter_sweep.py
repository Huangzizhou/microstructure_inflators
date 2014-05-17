#!/usr/bin/env python

import argparse
import json
import os.path
import numpy as np
from WireNetwork import WireNetwork

def load(wire_file):
    wire_network = WireNetwork();
    wire_network.load_from_file(wire_file);
    return wire_network;

def load_orbits(wire_file):
    basename, ext = os.path.splitext(wire_file);
    orbit_file = basename + ".orbit";
    with open(orbit_file, 'r') as fin:
        orbit_config = json.load(fin);
        return orbit_config["vertex_orbits"],\
                orbit_config["edge_orbits"];

def get_thickness_orbits(thickness_config, wire_network, vertex_orbits,
        edge_orbits):
    orbits = vertex_orbits if not thickness_config["per_edge"]\
            else edge_orbits;
    orbits = np.array(orbits);

    which_orbits = thickness_config["orbits"];
    if which_orbits == "all":
        pass;
    elif which_orbits == "none":
        orbits = [];
    else:
        orbits = orbits[which_orbits];
    return orbits;

def generate_default_thickness(thickness_config, wire_network):
    num_thicknesses = wire_network.num_vertices \
            if not thickness_config["per_edge"]\
            else wire_network.num_edges;
    thickness = np.ones(num_thicknesses) * thickness_config["default"];
    return thickness;

def param_to_thickness(param, wire_network, orbits, thickness_config):
    thickness = generate_default_thickness(thickness_config, wire_network);
    for value, orbit in zip(param, orbits):
        thickness[orbit] = value;
    return thickness;

def sweep_thickness(thickness_config, wire_network,
        vertex_orbits, edge_orbits):
    """
    return list of all possible thickness combinations + the parameter used.
    """
    orbits = get_thickness_orbits(thickness_config, wire_network,
            vertex_orbits, edge_orbits);
    num_dof = len(orbits);
    if num_dof != 0:
        thickness_range = thickness_config["thickness_range"];
        thickness_samples = thickness_config["num_samples"];
        values = np.linspace(thickness_range[0], thickness_range[1],
                thickness_samples);
        grid = np.meshgrid(*([values] * num_dof));
        parameters = np.vstack(map(np.ravel, grid)).T;

        thicknesses = [\
                param_to_thickness(param, wire_network, orbits, thickness_config) \
                for param in parameters];
    else:
        thicknesses = np.array([generate_default_thickness(thickness_config,
            wire_network)]);
        parameters = np.array([[]]);
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

def get_offset_orbits(thickness_config, wire_network, vertex_orbits):
    orbits = np.array(vertex_orbits);

    which_orbits = thickness_config["orbits"];
    if which_orbits == "all":
        pass;
    elif which_orbits == "none":
        orbits = [];
    else:
        orbits = orbits[which_orbits];
    return orbits;

def param_to_offset(param, wire_network, orbits):
    offset = np.zeros((wire_network.num_vertices, wire_network.dim));
    dof_count = 0;
    for orbit in orbits:
        orbit_num_dof = orbit.num_dof;
        orbit_dof_offset = orbit.get_vertex_offsets(
                param[dof_count:dof_count+orbit_num_dof]);
        offset[orbit.orbit] = orbit_dof_offset;
        dof_count += orbit_num_dof;
    return offset;

def sweep_offset(offset_config, wire_network, vertex_orbits):
    """
    return the list of all possible offset combinations + the parameter used.
    """
    orbits = get_offset_orbits(offset_config, wire_network, vertex_orbits);
    orbits = [VertexOrbit(wire_network, orbit) for orbit in orbits];
    num_orbits = len(orbits)
    num_dof = np.sum([orbit.num_dof for orbit in orbits]);

    percentage_range = offset_config["offset_percentage"];
    num_samples = offset_config["num_samples"];
    percentages = np.linspace(percentage_range[0], percentage_range[1], num_samples);
    offset_percent_space = [percentages] * num_dof;
    grid = np.meshgrid(*offset_percent_space);
    parameters = np.vstack(map(np.ravel, grid)).T;

    offsets = [param_to_offset(param, wire_network, orbits) for param in parameters];

    return offsets, parameters;

def save_setting(output_file, sweep_config, thickness_settings, offset_settings):
    thicknesses, thickness_parameters = thickness_settings;
    offsets, offset_parameters = offset_settings;

    if sweep_config["thickness"]["per_edge"]:
        thickness_key = "edge_thickness";
    else:
        thickness_key = "vertex_thickness";

    basename, ext = os.path.splitext(output_file);

    num_thicknesses = len(thicknesses);
    num_offsets = len(offsets);

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

def load_config(config_file):
    config_path = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);
        if not os.path.isabs(config["wire_file"]):
            config["wire_file"] = os.path.join(config_path, config["wire_file"]);
        return config;

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate parameter file");
    parser.add_argument("--output", "-o", help="output config files",
            required=True);
    parser.add_argument("sweep_file",
            help="configuration file of the sweep");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    sweep_config = load_config(args.sweep_file);

    wire_file = sweep_config["wire_file"];
    wire_network = load(wire_file);
    vertex_orbits, edge_orbits = load_orbits(wire_file);

    thickness_settings = sweep_thickness(sweep_config["thickness"],
            wire_network, vertex_orbits, edge_orbits);
    offset_settings = sweep_offset(sweep_config["offset"],
            wire_network, vertex_orbits);

    save_setting(args.output, sweep_config, thickness_settings, offset_settings);

if __name__ == "__main__":
    main();


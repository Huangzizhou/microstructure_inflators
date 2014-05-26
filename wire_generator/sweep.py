#!/usr/bin/env python

import argparse
import itertools
import json
import numpy as np
import os.path

from WireNetwork import WireNetwork

def load_config(config_file):
    config_path = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);
        if not os.path.isabs(config["wire_file"]):
            config["wire_file"] = os.path.join(config_path, config["wire_file"]);
        if not os.path.isabs(config["orbit_file"]):
            config["orbit_file"] = os.path.join(config_path, config["orbit_file"]);
        return config;

def load(wire_file):
    wire_network = WireNetwork();
    wire_network.load_from_file(wire_file);
    return wire_network;

def load_orbits(orbit_file):
    with open(orbit_file, 'r') as fin:
        orbit_config = json.load(fin);
        return np.array(orbit_config["vertex_orbits"]),\
                np.array(orbit_config["edge_orbits"]);

def sweep_thickness(thickness_config, wire_network):
    num_dof = len(thickness_config["orbits"]);
    if num_dof != 0:
        thickness_range = thickness_config["thickness_range"];
        thickness_samples = thickness_config["num_samples"];
        values = np.linspace(thickness_range[0], thickness_range[1],
                thickness_samples);
        grid = np.meshgrid(*([values] * num_dof));
        thicknesses = np.vstack(map(np.ravel, grid)).T;
    else:
        thicknesses = np.array([]);
    parameters = thicknesses;
    return thicknesses, parameters;

def compute_orbit_mask(wire_network, orbits):
    tol = 1e-6;
    mask = [];
    bbox_min, bbox_max = wire_network.bbox;
    for orbit in orbits:
        orbit_mask = np.array([False, False, False]);
        vertices = wire_network.vertices[orbit];
        orbit_bbox_min = np.amin(vertices, axis=0);
        orbit_bbox_max = np.amax(vertices, axis=0);
        orbit_bbox_size = orbit_bbox_max - orbit_bbox_min;

        zero_dim = np.absolute(orbit_bbox_size) < tol;
        on_min_boundary = np.absolute(orbit_bbox_min - bbox_min) < tol;
        on_max_boundary = np.absolute(orbit_bbox_max - bbox_max) < tol;

        orbit_mask = np.logical_or(orbit_mask, zero_dim);
        orbit_mask = np.logical_or(orbit_mask, on_min_boundary);
        orbit_mask = np.logical_or(orbit_mask, on_max_boundary);

        mask.append(orbit_mask);
    return np.array(mask);

def sweep_vertex_offset(offset_config, wire_network, vertex_orbits):
    orbits_indices = offset_config["orbits"];
    orbits_mask = compute_orbit_mask(wire_network, vertex_orbits[orbits_indices]);
    num_dof = np.count_nonzero(np.logical_not(orbits_mask));
    if num_dof != 0:
        offset_range = offset_config["offset_percentage"];
        offset_samples = offset_config["num_samples"];
        offset_default = offset_config["default"];
        values = np.linspace(offset_range[0], offset_range[1], offset_samples);
        grid = np.meshgrid(*([values] * num_dof));
        parameters = np.vstack(map(np.ravel, grid)).T;

        offsets = [];
        for param in parameters:
            offset = np.zeros_like(orbits_mask, dtype=float);
            offset[orbits_mask] = offset_default;
            offset[np.logical_not(orbits_mask)] = param;
            offsets.append(offset);
    else:
        offsets = np.array([]);
        parameters = offsets;
    return offsets, parameters;

def save_setting(output_file, sweep_config, thickness_settings, offset_settings):
    basename, ext = os.path.splitext(output_file);
    basedir = os.path.dirname(basename);
    orbit_file = sweep_config["orbit_file"];
    orbit_file = os.path.relpath(orbit_file, basedir);

    thicknesses, thickness_parameters = thickness_settings;
    vertex_offsets, offset_parameters = offset_settings;

    thickness_configs = [];
    if len(thicknesses) > 0:
        thickness_config = sweep_config["thickness"];
        thickness_type = "edge_orbit" if thickness_config["per_edge"]\
                else "vertex_orbit";
        for thickness in thicknesses:
            config = {
                    "type": thickness_type,
                    "orbit_file": orbit_file,
                    "effective_orbits": thickness_config["orbits"],
                    "thickness": thickness.tolist(),
                    "default": thickness_config["default"] };
            thickness_configs.append(config);

    offset_configs = [];
    if len(vertex_offsets) > 0:
        offset_config = sweep_config["offset"];
        for offset in vertex_offsets:
            config = {
                    "type": "vertex_orbit",
                    "orbit_file": orbit_file,
                    "effective_orbits": offset_config["orbits"],
                    "offset_percentages": offset.tolist() };
            offset_configs.append(config);

    for i,j in itertools.product(
            range(len(thickness_configs)),
            range(len(offset_configs))):
        thickness = thickness_configs[i];
        thickness_param = thickness_parameters[i].tolist();
        offset = offset_configs[j];
        offset_param = offset_parameters[j].tolist();
        contents = {
                "thickness": thickness,
                "thickness_parameter": thickness_param,
                "vertex_offset": offset,
                "offset_parameter": offset_param
                };
        filename = "{}_thickness_{}_offset_{}{}".format(basename, i, j, ext);
        with open(filename, 'w') as fout:
            json.dump(contents, fout, indent=4);

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
    orbit_file = sweep_config["orbit_file"];
    wire_network = load(wire_file);
    vertex_orbits, edge_orbits = load_orbits(orbit_file);

    thickness_setting = sweep_thickness(sweep_config["thickness"], wire_network);
    offset_setting = sweep_vertex_offset(sweep_config["offset"], wire_network,
            vertex_orbits);

    save_setting(args.output, sweep_config, thickness_setting, offset_setting);

if __name__ == "__main__":
    main();


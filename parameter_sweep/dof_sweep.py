#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path
from subprocess import check_call

import microstructures_setting
from wire_generator.core.WireNetwork import WireNetwork
from wire_generator.utils.find_file import find_file
from wire_generator.parameter.PyParameters import PyParameters
import PyWires

def load_sweep_config(config_file):
    """ syntax:
    {
        "dof_type": "isotropic" | "orthotropic",
        "wires": "path_to_wires.txt",
        "thickness_type": "vertex" | "edge",
        "cell_size": #,
        "thickness_range": [#, #],
        "thickness_samples": #,
        "offset_range": [#, #],
        "offset_samples": #
    }
    """
    with open(config_file, 'r') as fin:
        config = json.load(fin);
    root_dir = os.path.dirname(config_file);
    if not os.path.isabs(config["wires"]):
        wires_file = config["wires"];
        wires_file = os.path.join(root_dir, wires_file);
        wires_file = os.path.abspath(wires_file);
        config["wires"] = str(wires_file);
    return config;

def load_wires(wires_file):
    root_dir = os.path.dirname(wires_file);
    wires = [];
    with open(wires_file, 'r') as fin:
        for filename in fin:
            filename = filename.strip();
            if not os.path.isabs(filename):
                filename = find_file(filename, root_dir);
            assert(os.path.exists(filename));
            wire_network = WireNetwork();
            wire_network.load_from_file(filename);
            wire_network.filename = filename;
            wires.append(wire_network);
    return wires;

def get_pattern_name(wire_network):
    wire_file = wire_network.filename;
    basename = os.path.basename(wire_file);
    name, ext = os.path.splitext(basename);
    return name

def copy_wire_to_out_dir(wire_network, out_dir):
    wire_file = wire_network.filename;
    local_wire_file = os.path.join(out_dir, os.path.basename(wire_file));
    command = "cp {} {}".format(wire_file, local_wire_file);
    check_call(command.split());
    return local_wire_file;

def generate_config_file(config_file_name, dof_file_name, wire_network,
        cell_size, out_dir):
    local_wire_file = copy_wire_to_out_dir(wire_network, out_dir);
    bbox_min = -np.ones(wire_network.dim) * 0.5 * cell_size;
    bbox_max =  np.ones(wire_network.dim) * 0.5 * cell_size;
    repeats = np.ones(wire_network.dim);
    config = {
            "periodic": True,
            "thickness": 0.5,
            "wire_network": os.path.relpath(local_wire_file, out_dir),
            "dof_file": os.path.relpath(dof_file_name, out_dir),
            "bbox_min": bbox_min.tolist(),
            "bbox_max": bbox_max.tolist(),
            "repeats": repeats.tolist()
            };
    with open(config_file_name, 'w') as fout:
        json.dump(config, fout, indent=4);

def generate_sweeps(wire_network, config):
    if config["thickness_type"] == "vertex":
        thickness_type = PyWires.VERTEX;
    else:
        thickness_type = PyWires.EDGE;

    parameters = PyParameters(wire_network, 0.5);
    if config["dof_type"] == "isotropic":
        parameters.load_default_isotropic_parameters(thickness_type);
    else:
        parameters.load_default_orthotropic_parameters(thickness_type);
    num_thickness_dofs = parameters.num_thickness_dofs;
    num_offset_dofs = parameters.num_offset_dofs;

    thickness_samples = np.linspace(
            config["thickness_range"][0],
            config["thickness_range"][1],
            config["thickness_samples"]);
    offset_samples = np.linspace(
            config["offset_range"][0],
            config["offset_range"][1],
            config["offset_samples"]);

    dof_samples = [thickness_samples for i in range(num_thickness_dofs) ] +\
            [offset_samples for i in range(num_offset_dofs)] ;
    dofs = np.meshgrid(*dof_samples);
    dofs = np.array([dof.ravel(order="C") for dof in dofs]).T;

    pattern_name = get_pattern_name(wire_network);
    print("{}: {:3} dofs, {:3} thickness, {:3} offset, {:6} combinations".format(
        pattern_name, parameters.num_dofs,
        num_thickness_dofs, num_offset_dofs,
        len(dofs)));
    if config["dry_run"]: return len(dofs);

    config_files = [];
    for i,dof_sample in enumerate(dofs):
        basename = "{}_sample_{:06}".format(pattern_name, i);
        basename = os.path.join(config["out_dir"], basename);
        dof_file_name = "{}.dof".format(basename);
        parameters.dofs = dof_sample;
        parameters.save_dof_file(dof_file_name);

        config_file_name = "{}.config".format(basename);
        generate_config_file(
                config_file_name,
                dof_file_name,
                wire_network,
                config["cell_size"],
                config["out_dir"]);
        config_files.append(config_file_name);

    return config_files;

def output_job_file(output_dir, all_config_files):
    job_file = os.path.join(output_dir, "jobs.txt");
    with open(job_file, 'w') as fout:
        for name in all_config_files:
            name = os.path.relpath(name, output_dir);
            fout.write("{}\n".format(name));

def parse_args():
    parser = argparse.ArgumentParser(
            description="Generate dofs files for isotropic sweep");
    parser.add_argument("--out-dir", "-o", help="output directory",
            required=True);
    parser.add_argument("--dry-run", help="run without generating dof files",
            action="store_true");
    parser.add_argument("sweep_config_file");
    return parser.parse_args();

def main():
    args = parse_args();
    config = load_sweep_config(args.sweep_config_file);
    config["out_dir"] = args.out_dir;
    config["dry_run"] = args.dry_run;
    wires = load_wires(config["wires"]);
    print("{} patterns loaded".format(len(wires)));

    all_config_files = [];
    for wire in wires:
        config_files = generate_sweeps(wire, config);
        all_config_files += config_files;
    print("{} total dof values sampled".format(len(all_config_files)));
    output_job_file(args.out_dir, all_config_files);

if __name__ == "__main__":
    main();

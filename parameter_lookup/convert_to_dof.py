#!/usr/bin/env python

import argparse
import json
import numpy as np
from numpy.linalg import norm
import os.path

import microstructures_setting
from wire_generator.core.WireNetwork import WireNetwork
from wire_generator.parameter.PyParameters import PyParameters
from wire_generator.inflator.PyWiresTiler import PyWiresTiler
import PyWires

def load_wire_network(filename):
    assert(os.path.exists(filename));
    wire_network = WireNetwork();
    wire_network.load_from_file(filename);
    wire_network.filename = filename;
    return wire_network;

def extract_thickness_dofs_from_modifier_file(wire_network, modifier_file):
    assert(os.path.exists(modifier_file));
    parameters = PyParameters(wire_network, 0.5);
    parameters.load_modifier_file(modifier_file);
    num_thickness_dofs = parameters.num_thickness_dofs;
    return parameters.dofs[:num_thickness_dofs];

def extract_offset_dofs_from_wires(parameters, source_wires, target_wires):
    assert(source_wires.num_vertices == target_wires.num_vertices);

    num_dofs = parameters.num_dofs;
    num_thickness_dofs = parameters.num_thickness_dofs;
    num_offset_dofs = parameters.num_offset_dofs;

    offset = target_wires.vertices - source_wires.vertices;

    offset_param = np.zeros(parameters.num_offset_dofs);
    for i in range(num_thickness_dofs, num_dofs):
        grad = parameters.raw_parameters.compute_wire_gradient(i);
        if np.all(grad == 0.0): continue;
        non_zero_entries = grad != 0.0;
        param_value = np.divide(
                offset[non_zero_entries], grad[non_zero_entries]);
        ave_value = np.mean(param_value);
        assert(np.amax(np.absolute(param_value - ave_value)) < 1e-12);
        offset_param[i - num_thickness_dofs] = ave_value;
    return offset_param;

def check_dof(source_wires, target_wires, parameters):
    tiler = PyWiresTiler();
    tiler.set_single_cell_from_wire_network(source_wires);
    bbox_min, bbox_max = source_wires.bbox;
    repeats = np.ones(3, dtype=int);
    tiler.tile(bbox_min, bbox_max, repeats, parameters);
    out_wires = tiler.wire_network;
    out_wires.translate(source_wires.bbox_center - out_wires.bbox_center);
    assert(norm(out_wires.vertices - target_wires.vertices) < 1e-6);

def compute_dofs(source_wires, target_wires, modifier_file):
    parameters = PyParameters(source_wires, 0.5);
    parameters.load_default_isotropic_parameters(PyWires.EDGE);
    dofs = parameters.dofs;

    if modifier_file is not None:
        thickness_dofs = extract_thickness_dofs_from_modifier_file(
                source_wires, modifier_file);
        assert(len(thickness_dofs) == parameters.num_thickness_dofs);
        dofs[:len(thickness_dofs)] = thickness_dofs;

    offset_dofs = extract_offset_dofs_from_wires(
            parameters, source_wires, target_wires);
    dofs[len(thickness_dofs):] = offset_dofs;
    parameters.dofs = dofs;
    check_dof(source_wires, target_wires, parameters);
    return parameters;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Convert offsetted wire_network and modifier into dof file");
    parser.add_argument("--modifier-file",
            help="Modifier file for offset wire network", default=None);
    parser.add_argument("reference_wire_file", help="reference wire network");
    parser.add_argument("offset_wire_file", help="offset wire network");
    parser.add_argument("output_file", help="output file");
    return parser.parse_args();

def main():
    args = parse_args();

    reference_wires = load_wire_network(args.reference_wire_file);
    offset_wires = load_wire_network(args.offset_wire_file);

    parameters = compute_dofs(reference_wires, offset_wires, args.modifier_file);
    parameters.save_dof_file(args.output_file);

if __name__ == "__main__":
    main();

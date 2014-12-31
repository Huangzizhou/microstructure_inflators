#!/usr/bin/env python

import argparse
import os.path
import numpy as np

import microstructures_setting
from wire_generator.utils.find_file import find_file
from wire_generator.tile_only import generate_box_guide_mesh
from wire_generator.tile import load_mesh, save_mesh, load_wire
from wire_generator.parameter.PyParameters import PyParameters

def load_wires(wire_list_file):
    root_dir = os.path.dirname(wire_list_file);

    with open(wire_list_file, 'r') as fin:
        wire_files = [str(filename.strip()) for filename in fin];

    wires = [];
    for name in wire_files:
        if not os.path.isabs(name):
            name = find_file(name, root_dir);
        wires.append(load_wire(name));
    return wires;

def compute_default_dofs(wires, thickness):
    all_dofs = [];
    for wire_network in wires:
        params = PyParameters(wire_network, thickness);
        params.load_default_isotropic_parameters();
        all_dofs.append(params.dofs);
    return all_dofs;

def generate_guide_mesh(dim, num_cells, cell_size):
    side_len = num_cells * cell_size;
    bbox_min = np.ones(dim) * (-0.5 * side_len);
    bbox_max = np.ones(dim) * (0.5 * side_len);
    num_cells = np.ones(dim, dtype=int) * num_cells;

    mesh = generate_box_guide_mesh(dim, bbox_min, bbox_max, num_cells);
    mesh.add_attribute("voxel_index");
    return mesh;

def assign_dofs(mesh, all_dofs):
    num_patterns = len(all_dofs);
    if mesh.get_dim() == 3:
        num_cells = mesh.get_num_voxels();
    else:
        num_cells = mesh.get_num_faces();
    pattern_indices = np.arange(num_cells, dtype=int);
    pattern_indices %= num_patterns;
    mesh.add_attribute("pattern_id");
    mesh.set_attribute("pattern_id", pattern_indices.astype(float));

    max_num_dofs = np.amax([len(dofs) for dofs in all_dofs]);
    dofs_table = np.zeros((num_cells, max_num_dofs));
    for i,dofs in enumerate(all_dofs):
        dofs_table[pattern_indices == i, :len(dofs)] = dofs;

    for i in range(max_num_dofs):
        attr_name = "dof_{}".format(i);
        attr_val = dofs_table[:, i].ravel();
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, attr_val);

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate hex mesh with checker board pattern");
    parser.add_argument("--num-cells", type=int,
            help="number of cells in each dimension", default=3);
    parser.add_argument("--cell-size", help="cell size in mm",
            type=int, default=5);
    parser.add_argument("--thickness", help="thickness value in mm",
            type=float, default=0.5);
    parser.add_argument("wire_list_file", help="wire list file");
    parser.add_argument("output_file", help="output mesh file");
    return parser.parse_args();

def main():
    args = parse_args();
    all_wires = load_wires(args.wire_list_file);
    all_dofs = compute_default_dofs(all_wires, args.thickness);
    mesh = generate_guide_mesh(
            all_wires[0].dim, args.num_cells, args.cell_size);
    assign_dofs(mesh, all_dofs);

    save_mesh(mesh, args.output_file);

if __name__ == "__main__":
    main();


#!/usr/bin/env python

import argparse
import json
import os.path
import numpy as np
from subprocess import check_call

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
        network = load_wire(name);
        network.filename = name;
        wires.append(network);
    return wire_files, wires;

def compute_default_dofs(wires, thickness, output_dir):
    all_dofs = [];
    dof_files = [];
    for wire_network in wires:
        basename = os.path.basename(wire_network.filename);
        name, ext = os.path.splitext(basename);
        dof_file = os.path.join(output_dir, "{}.dof".format(name));
        params = PyParameters(wire_network, thickness);
        params.load_default_isotropic_parameters();
        params.raw_parameters.save_dofs(dof_file);
        all_dofs.append(params.dofs);
        dof_files.append(dof_file);
    return all_dofs, dof_files;

def generate_config_files(dim, wire_files, dof_files, cell_size):
    bbox_min = np.zeros(dim).tolist();
    bbox_max = (np.ones(dim) * cell_size).tolist();
    repeats = np.ones(dim).tolist();
    config_files = [];
    for wire_file, dof_file in zip(wire_files, dof_files):
        name, ext = os.path.splitext(dof_file);
        config_file = "{}.config".format(name);
        config = {
                "thickness": 0.5,
                "periodic": True,
                "wire_network": wire_file,
                "dof_file": dof_file,
                "bbox_min": bbox_min,
                "bbox_max": bbox_max,
                "repeats": repeats
                }
        with open(config_file, 'w') as fout:
            json.dump(config, fout, indent=4);

        config_files.append(config_file);
    return config_files;

def generate_guide_mesh(dim, num_cells, cell_size):
    side_len = num_cells * cell_size;
    bbox_min = np.ones(dim) * (-0.5 * side_len);
    bbox_max = np.ones(dim) * (0.5 * side_len);
    num_cells = np.ones(dim, dtype=int) * num_cells;

    mesh = generate_box_guide_mesh(dim, bbox_min, bbox_max, num_cells);
    return mesh;

def get_num_cells(mesh):
    if mesh.get_dim() == 3:
        num_cells = mesh.get_num_voxels();
    else:
        num_cells = mesh.get_num_faces();
    return num_cells;

def assign_dofs(mesh, all_dofs):
    num_patterns = len(all_dofs);
    num_cells = get_num_cells(mesh);
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

def assign_materials(mesh, materials):
    pattern_id = mesh.get_attribute("pattern_id").ravel().astype(int);
    num_cells = get_num_cells(mesh);

    if mesh.get_dim() == 3:
        young = np.zeros((num_cells, 3));
        shear = np.zeros((num_cells, 3));
        poisson = np.zeros((num_cells, 6));
        young_names = ["young_x", "young_y", "young_z"];
        shear_names = ["shear_yz", "shear_zx", "shear_xy"];
        poisson_names = [
                "poisson_yz", "poisson_zy",
                "poisson_zx", "poisson_xz",
                "poisson_xy", "poisson_yx", ];
    else:
        young = np.zeros((num_cells, 2));
        shear = np.zeros((num_cells, 1));
        poisson = np.zeros((num_cells, 2));
        young_names = ["young_x", "young_y"];
        shear_names = ["shear_xy"];
        poisson_names = ["poisson_xy", "poisson_yx", ];

    for i,pattern_index in enumerate(pattern_id):
        young[i] = materials[pattern_index]["young"];
        shear[i] = materials[pattern_index]["shear"];
        poisson[i] = materials[pattern_index]["poisson"];

    for i,name in enumerate(young_names):
        mesh.add_attribute(name);
        mesh.set_attribute(name, young[:,i].ravel());
    for i,name in enumerate(shear_names):
        mesh.add_attribute(name);
        mesh.set_attribute(name, shear[:,i].ravel());
    for i,name in enumerate(poisson_names):
        mesh.add_attribute(name);
        mesh.set_attribute(name, poisson[:,i].ravel());

    mixed_material = {
            "type": "element_wise_orthotropic_material",
            "young": young_names,
            "shear": shear_names,
            "poisson": poisson_names,
            "density": 1.0
            };

    with open("mixed.material", 'w') as fout:
        json.dump(mixed_material, fout, indent=4);

def run_homogenization(config_files):
    exe_name = os.path.join(microstructures_setting.MICROSTRUCTURES_PATH,
            "wire_generator/tile_and_fit_periodic.py");
    material_file = os.path.join(microstructures_setting.MICROSTRUCTURES_PATH,
            "wire_generator/examples/Projet7000.material");
    materials = [];
    for config_file in config_files:
        name, ext = os.path.splitext(config_file);
        output_file = "{}.msh".format(name);
        command = "{} --material {} {} {}".format(exe_name, material_file,
            config_file, output_file);
        check_call(command.split());

        out_material_file = "{}.material".format(name);
        with open(out_material_file, 'r') as fin:
            material = json.load(fin);
            materials.append(material);
    return materials;

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
    output_dir = os.path.dirname(args.output_file);

    wire_files, all_wires = load_wires(args.wire_list_file);
    dim = all_wires[0].dim;
    all_dofs, dof_files = compute_default_dofs(all_wires, args.thickness, output_dir);
    config_files = generate_config_files(dim, wire_files, dof_files, args.cell_size);
    mesh = generate_guide_mesh(
            dim, args.num_cells, args.cell_size);
    assign_dofs(mesh, all_dofs);

    homogenized_materials = run_homogenization(config_files);
    assign_materials(mesh, homogenized_materials);

    save_mesh(mesh, args.output_file);

if __name__ == "__main__":
    main();


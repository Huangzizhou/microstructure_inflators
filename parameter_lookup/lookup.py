#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

import PyMeshSetting
import PyMesh
from IsotropicMaterial import IsotropicMaterial
from PatternParameterTable import PatternParameterTable

def load_json(filename):
    with open(filename, 'r') as fin:
        config = json.load(fin);
        return config;

def load_mesh(mesh_file):
    factory = PyMesh.MeshFactory();
    factory.load_file(mesh_file);
    factory.drop_zero_dim();
    return factory.create();

def save_mesh(mesh_file, mesh, *attributes):
    writer = PyMesh.MeshWriter.create_writer(mesh_file);
    for attr in attributes:
        assert(mesh.has_attribute(attr));
        writer.with_attribute(attr);
    writer.write_mesh(mesh);

def extract_material_properties(mesh):
    young = mesh.get_attribute("young").ravel();
    poisson = mesh.get_attribute("poisson").ravel();
    return young, poisson;

def output_config_file(guide_mesh_file, index_dir):
    name, ext = os.path.splitext(guide_mesh_file);
    sweep_file = os.path.join(index_dir, "../wires.sweep");
    root_dir = os.path.dirname(sweep_file);
    assert(os.path.exists(sweep_file));
    sweep_config = load_json(sweep_file);
    if not os.path.isabs(sweep_config["wires"]):
        sweep_config["wires"] = os.path.join(
                root_dir, sweep_config["wires"]);
    config = {
            "wire_list_file": sweep_config["wires"],
            "guide_mesh": guide_mesh_file,
            "subdiv": 0,
            "dof_type": sweep_config["dof_type"],
            "thickness_type": sweep_config["thickness_type"],
            }
    config_file = "{}.config".format(name);
    with open(config_file, 'w') as fout:
        json.dump(config, fout, indent=4);

def parse_args():
    parser = argparse.ArgumentParser(
            description="lookup pattern parameter from material properties");
    parser.add_argument("--metric", help="metric of proximity",
            default="compliance", choices=("compliance", "elasticity"));
    parser.add_argument("--index-dir", help="index directory");
    parser.add_argument("--rehomogenize", action="store_true",
            help="rerun homogenization on fitted pattern parameters");
    parser.add_argument("--material", default=None,
            help="material filed used for homogenization");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    if args.rehomogenize and args.material is None:
        raise RuntimeError("Please specify a material file to use in homogenization");

    mesh = load_mesh(args.input_mesh);

    young, poisson = extract_material_properties(mesh);

    materials = [IsotropicMaterial(mesh.get_dim(), E, nu)
            for E, nu in zip(young, poisson)];

    param_table = PatternParameterTable(
            args.index_dir, args.metric, args.material);
    header = param_table.header;

    param_values, material_header, material_values =\
            param_table.lookup_and_interpolate(materials,
                    rehomogenize=args.rehomogenize,
                    num_candidates=1);
    for i,attr_name in enumerate(header):
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, param_values[:,i]);

    for i,attr_name in enumerate(material_header):
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, material_values[:,i]);

    header += material_header;
    header.append("young");
    header.append("poisson");

    save_mesh(args.output_mesh, mesh, *header);
    output_config_file(args.output_mesh, args.index_dir);

if __name__ == "__main__":
    main();


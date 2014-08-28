#!/usr/bin/env python

import argparse
import numpy as np

import PyMeshSetting
import PyMesh
from IsotropicMaterial import IsotropicMaterial
from PatternParameterTable import PatternParameterTable

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

def parse_args():
    parser = argparse.ArgumentParser(
            description="lookup pattern parameter from material properties");
    parser.add_argument("--index-dir", help="index directory");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);

    young, poisson = extract_material_properties(mesh);

    materials = [IsotropicMaterial(mesh.get_dim(), E, nu)
            for E, nu in zip(young, poisson)];

    param_table = PatternParameterTable(args.index_dir);
    header = param_table.header;

    param_values, dist = param_table.lookup(materials);
    for i,attr_name in enumerate(header):
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, param_values[:,i]);

    save_mesh(args.output_mesh, mesh, *header);
    
if __name__ == "__main__":
    main();


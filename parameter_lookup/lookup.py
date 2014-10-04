#!/usr/bin/env python

import argparse
import numpy as np
import os.path

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
    parser.add_argument("--metric", help="metric of proximity",
            default="compliance", choices=("compliance", "elasticity"));
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

    param_table = PatternParameterTable(args.index_dir, args.metric);
    header = param_table.header;

    #param_values, fitted_young, fitted_poisson, fitted_shear, dist =\
    #        param_table.lookup(materials);
    param_values, fitted_young, fitted_poisson, fitted_shear, dist =\
            param_table.lookup_and_interpolate(materials);
    for i,attr_name in enumerate(header):
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, param_values[:,i]);

    header.append("young");
    header.append("poisson");
    mesh.add_attribute("fitted_young_x");
    mesh.set_attribute("fitted_young_x", fitted_young[:,0]);
    header.append("fitted_young_x");
    mesh.add_attribute("fitted_young_y");
    mesh.set_attribute("fitted_young_y", fitted_young[:,1]);
    header.append("fitted_young_y");
    if mesh.get_dim() == 2:
        mesh.add_attribute("fitted_poisson_xy");
        mesh.set_attribute("fitted_poisson_xy", fitted_poisson[:,0]);
        mesh.add_attribute("fitted_poisson_yx");
        mesh.set_attribute("fitted_poisson_yx", fitted_poisson[:,1]);
        header.append("fitted_poisson_xy");
        header.append("fitted_poisson_yx");
        mesh.add_attribute("shear_xy");
        mesh.set_attribute("shear_xy", fitted_shear[:,0]);
        header.append("shear_xy");
    elif mesh.get_dim() == 3:
        mesh.add_attribute("fitted_young_z");
        mesh.set_attribute("fitted_young_z", fitted_young[:,2]);
        header.append("fitted_young_z");

        mesh.add_attribute("fitted_poisson_yx");
        mesh.set_attribute("fitted_poisson_yz", fitted_poisson[:,0]);
        mesh.add_attribute("fitted_poisson_zy");
        mesh.set_attribute("fitted_poisson_zy", fitted_poisson[:,1]);
        mesh.add_attribute("fitted_poisson_zx");
        mesh.set_attribute("fitted_poisson_zx", fitted_poisson[:,2]);
        mesh.add_attribute("fitted_poisson_xz");
        mesh.set_attribute("fitted_poisson_xz", fitted_poisson[:,3]);
        mesh.add_attribute("fitted_poisson_xy");
        mesh.set_attribute("fitted_poisson_xy", fitted_poisson[:,4]);
        mesh.add_attribute("fitted_poisson_yx");
        mesh.set_attribute("fitted_poisson_yx", fitted_poisson[:,5]);
        header = header + [
                "fitted_poisson_yz",
                "fitted_poisson_zy",
                "fitted_poisson_zx",
                "fitted_poisson_xz",
                "fitted_poisson_yx",
                "fitted_poisson_xy" ];
        mesh.add_attribute("shear_yz");
        mesh.set_attribute("shear_yz", fitted_shear[:,0]);
        mesh.add_attribute("shear_zx");
        mesh.set_attribute("shear_zx", fitted_shear[:,1]);
        mesh.add_attribute("shear_xy");
        mesh.set_attribute("shear_xy", fitted_shear[:,2]);
        header.append("shear_yz");
        header.append("shear_zx");
        header.append("shear_xy");

    save_mesh(args.output_mesh, mesh, *header);

if __name__ == "__main__":
    main();


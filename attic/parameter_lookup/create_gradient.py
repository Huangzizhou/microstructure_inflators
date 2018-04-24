#!/usr/bin/env python

import argparse
import numpy as np
import os.path

import PyMeshSetting
import PyMesh

from IsotropicMaterial import IsotropicMaterial
from PatternParameterTable import PatternParameterTable

AXIS = ["X", "Y", "Z"];

def form_mesh(vertices, faces, voxels=np.array([])):
    factory = PyMesh.MeshFactory();
    factory.load_data(
            vertices.ravel(order="C"),
            faces.ravel(order="C"),
            voxels.ravel(order="C"), 2, 4, 8);
    return factory.create();

def save_mesh(mesh_file, mesh, *attributes):
    writer = PyMesh.MeshWriter.create_writer(mesh_file);
    for attr in attributes:
        assert(mesh.has_attribute(attr));
        writer.with_attribute(attr);
    writer.write_mesh(mesh);

def generate_quad_mesh(rows, cols, cell_size):
    grid_x = np.linspace(0.0, cols * cell_size, cols+1);
    grid_y = np.linspace(0.0, rows * cell_size, rows+1);
    grid_x, grid_y = np.meshgrid(grid_x, grid_y);
    grid_x = grid_x.ravel();
    grid_y = grid_y.ravel();

    vertices = np.array([[x, y] for x,y in zip(grid_x, grid_y)]);
    step = cols+1;
    faces = np.array([
            [i*step+j, (i+1)*step+j, (i+1)*step+j+1, i*step+j+1]
            for i in range(rows)
            for j in range(cols)
            ], dtype=int);

    row_indices = np.repeat(np.arange(rows, dtype=int), cols);
    col_indices = np.tile(np.arange(cols, dtype=int), rows);

    mesh = form_mesh(vertices, faces);

    mesh.add_attribute("row_index");
    mesh.set_attribute("row_index", row_indices);
    mesh.add_attribute("col_index");
    mesh.set_attribute("col_index", col_indices);

    return mesh;

def get_young_range(index_dir):
    young_file = os.path.join(index_dir, "young.npy");
    young = np.load(young_file);
    min_young = np.amin(young);
    max_young = np.amax(young);
    young_range = max_young - min_young;
    #return min_young + 0.25 * young_range, max_young - 0.25 * young_range;
    return min_young, max_young - 0.1 * young_range;

def get_poisson_range(index_dir):
    poisson_file = os.path.join(index_dir, "poisson.npy");
    poisson = np.load(poisson_file);
    return np.amin(poisson), np.amax(poisson);

def create_material_gradient(mesh, index_dir):
    row_indices = mesh.get_attribute("row_index").ravel();
    col_indices = mesh.get_attribute("col_index").ravel();

    rows = np.amax(row_indices);
    cols = np.amax(col_indices);

    young_min, young_max = get_young_range(index_dir);
    poisson_min, poisson_max = get_poisson_range(index_dir);

    poisson_max = min(poisson_max, 0.8);
    poisson_min = max(poisson_min, 0.3);

    young_values = [];
    poisson_values = [];
    for i in range(mesh.get_num_faces()):
        row = row_indices[i];
        col = col_indices[i];

        #young = float(row) / rows * (young_max - young_min) + young_min;
        #young = 0.5 * (young_max + young_min);
        young = young_min + 0.5 * (young_max - young_min);
        #poisson = float(col) / cols * (poisson_max - poisson_min) + poisson_min;
        poisson = float(row) / rows * (poisson_max - poisson_min) + poisson_min;
        #poisson = 0.3;

        young_values.append(young);
        poisson_values.append(poisson);

    mesh.add_attribute("young");
    mesh.set_attribute("young", np.array(young_values));
    mesh.add_attribute("poisson");
    mesh.set_attribute("poisson", np.array(poisson_values));

def add_fitted_material_properties(mesh, fitted_young, fitted_poisson,
        fitted_shear):
    if mesh.get_dim() == 2:
        mesh.add_attribute("fitted_young_x");
        mesh.set_attribute("fitted_young_x", fitted_young[:,0]);
        mesh.add_attribute("fitted_young_y");
        mesh.set_attribute("fitted_young_y", fitted_young[:,1]);
        mesh.add_attribute("fitted_poisson_xy");
        mesh.set_attribute("fitted_poisson_xy", fitted_poisson[:,0]);
    elif mesh.get_dim() == 3:
        mesh.add_attribute("fitted_young_x");
        mesh.set_attribute("fitted_young_x", fitted_young[:,0]);
        mesh.add_attribute("fitted_young_y");
        mesh.set_attribute("fitted_young_y", fitted_young[:,1]);
        mesh.add_attribute("fitted_young_z");
        mesh.set_attribute("fitted_young_z", fitted_young[:,2]);

def lookup_pattern_parameter(mesh, index_dir):
    young = mesh.get_attribute("young");
    poisson = mesh.get_attribute("poisson");

    materials = [IsotropicMaterial(mesh.get_dim(), E, nu)
            for E, nu in zip(young, poisson)];

    param_table = PatternParameterTable(index_dir);
    header = param_table.header;

    #param_values, fitted_young, fitted_poisson, fitted_shear, dist =\
    #        param_table.lookup(materials);
    param_values, fitted_young, fitted_poisson, fitted_shear, dist =\
            param_table.lookup_and_interpolate(materials);
    for i,attr_name in enumerate(header):
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, param_values[:,i]);

    mesh.add_attribute("fit_error");
    mesh.set_attribute("fit_error", dist);

    add_fitted_material_properties(mesh, fitted_young, fitted_poisson,
            fitted_shear);


def parse_args():
    parser = argparse.ArgumentParser(
            description="Generating a regular grid of patterns based on changing material");
    parser.add_argument("--rows", type=int, help="number of rows", default=10);
    parser.add_argument("--cols", type=int, help="number of cols", default=10);
    parser.add_argument("--cell-size", type=float, help="cell size in mm", default=5);
    parser.add_argument("--index-dir", help="index directory");
    parser.add_argument("output_mesh", help="output mesh file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = generate_quad_mesh(args.rows, args.cols, args.cell_size);
    create_material_gradient(mesh, args.index_dir);
    lookup_pattern_parameter(mesh, args.index_dir);

    save_mesh(args.output_mesh, mesh, *mesh.get_attribute_names());

if __name__ == "__main__":
    main();


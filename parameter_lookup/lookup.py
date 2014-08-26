#!/usr/bin/env python

import argparse
import csv
import numpy as np
from numpy.linalg import norm
import os.path

import PyMeshSetting
import PyMesh

import pyflann

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

def load_index(index_dir, index_name):
    data = load_dataset(index_dir, index_name);
    table = pyflann.FLANN();
    table.load_index(os.path.join(index_dir, index_name + ".index"), data);
    return table;

def load_dataset(index_dir, dataset_name):
    data = np.load(os.path.join(index_dir, dataset_name + ".npy"));
    return data;

def load_data(csv_file):
    data = [];
    with open(csv_file, 'r') as fin:
        reader = csv.reader(fin);
        header = reader.next();
        for row in reader:
            data.append(row);
    return header, np.array(data);

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

    if mesh.get_dim() == 3:
        young = np.repeat(young, 3).reshape((-1, 3));
        poisson = np.repeat(poisson, 6).reshape((-1, 6));
    elif mesh.get_dim() == 2:
        young = np.repeat(young, 2).reshape((-1, 2));
        poisson = np.repeat(poisson, 2).reshape((-1, 2));

    young_table = load_index(args.index_dir, "young");
    poisson_table = load_index(args.index_dir, "poisson");
    poisson_data = load_dataset(args.index_dir, "poisson");
    header, pattern = load_data(os.path.join(args.index_dir, "pattern.csv"));

    young_index, young_dist = young_table.nn_index(young, 5);
    poisson_index, poisson_dist = poisson_table.nn_index(poisson, 5);

    index = [];
    for i, top_indices in enumerate(young_index):
        inferred_poisson = poisson_data[top_indices];
        poisson_dist = norm(inferred_poisson - poisson[i], axis=1);
        idx = np.argmin(poisson_dist);
        index.append(top_indices[idx]);

    #index = young_index[:, 0].ravel();
    for i,attr_name in enumerate(header):
        values = pattern[index, i].ravel();
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, values);

    save_mesh(args.output_mesh, mesh, *header);
    
if __name__ == "__main__":
    main();


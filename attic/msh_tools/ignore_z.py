#!/usr/bin/env python

import argparse
import numpy as np

from add_element_attribute import load_mesh, save_mesh

import PyMeshSetting
import PyMesh

def form_mesh(vertices, faces, voxels, dim, vertex_per_face, vertex_per_voxel):
    factory = PyMesh.MeshFactory();
    factory.load_data(vertices.ravel(), faces.ravel(), voxels.ravel(),
            dim, vertex_per_face, vertex_per_voxel);
    return factory.create();

def parse_args():
    parser = argparse.ArgumentParser(description="Discard z coordinates");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    return parser.parse_args();

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);

    if mesh.get_dim() == 2:
        print("mesh is already 2D.  Nothing to do...");
        return;

    vertices = mesh.get_vertices().reshape((-1, 3));
    faces = mesh.get_faces();
    voxels = mesh.get_voxels();

    vertices = vertices[:,0:2].ravel(order="C");
    mesh2 = form_mesh(vertices, faces, voxels,
            2, mesh.get_vertex_per_face(), mesh.get_vertex_per_voxel());

    attribute_names = mesh.get_attribute_names();
    for attr_name in attribute_names:
        mesh2.add_attribute(attr_name);
        mesh2.set_attribute(attr_name, mesh.get_attribute(attr_name).ravel());
    save_mesh(args.output_mesh, mesh2, *attribute_names);

if __name__ == "__main__":
    main();

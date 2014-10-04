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
    parser = argparse.ArgumentParser(description="scale mesh while keep all attirbutes");
    parser.add_argument("--scale", "-s", help="scale factor", type=float);
    parser.add_argument("--center", "-c",
            help="center scaled mesh at bbox center", action="store_true");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);

    vertices = mesh.get_vertices().reshape((mesh.get_num_vertices(), -1),
            order="C");
    vertices *= args.scale;

    if args.center:
        bbox_min = np.amin(vertices, axis=0);
        bbox_max = np.amax(vertices, axis=0);
        bbox_center = 0.5 * (bbox_min + bbox_max);
        vertices -= bbox_center;

    mesh2 = form_mesh(vertices, mesh.get_faces(), mesh.get_voxels(),
            mesh.get_dim(), mesh.get_vertex_per_face(),
            mesh.get_vertex_per_voxel());

    attribute_names = mesh.get_attribute_names();
    for attr_name in attribute_names:
        mesh2.add_attribute(attr_name);
        mesh2.set_attribute(attr_name, mesh.get_attribute(attr_name).ravel());

    save_mesh(args.output_mesh, mesh2, *attribute_names);

if __name__ == "__main__":
    main();

#!/usr/bin/env python

import argparse
import numpy as np
import LinearElasticitySettings
from BoxMeshGenerator import split_hex_into_tets, split_hex_into_tets_symmetrically
from BoxMeshGenerator import remove_duplicated_vertices, remove_isolated_vertices

import PyMesh

def merge_identical_vertices(vertices, voxels):
    vertices, voxels = remove_duplicated_vertices(vertices, voxels);
    vertices, voxels = remove_isolated_vertices(vertices, voxels);
    return vertices, voxels;

def load_mesh(filename):
    factory = PyMesh.MeshFactory();
    factory.load_file(filename);
    return factory.create();

def form_mesh(vertices, faces, voxels=np.array([])):
    dim = vertices.shape[1];
    factory = PyMesh.MeshFactory();
    if dim == 3:
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels.ravel(order="C"), 3, 3, 4);
    elif dim == 2:
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels.ravel(order="C"), 2, 3, 4);
    return factory.create();

def save_mesh(filename, mesh, *attributes):
    writer = PyMesh.MeshWriter.create_writer(filename);
    for attr in attributes:
        if not mesh.has_attribute(attr):
            raise KeyError("Attribute {} is not found in mesh".format(attr));
        writer.with_attribute(attr);
    writer.write_mesh(mesh);

def hex2tet(hex_file, tet_file, keep_symmetry):
    hex_mesh = load_mesh(hex_file);
    hex_vertices = hex_mesh.get_vertices().reshape(
            (hex_mesh.get_num_vertices(), -1));
    hexes = hex_mesh.get_voxels().reshape((-1, hex_mesh.get_vertex_per_voxel()));

    vertex_count = 0;
    vertices = [];
    voxels = [];

    for voxel in hexes:
        corners = hex_vertices[voxel];
        if keep_symmetry:
            vts, tets = split_hex_into_tets_symmetrically(corners);
        else:
            vts, tets = split_hex_into_tets(corners);
        vertices.append(vts);
        voxels.append(tets + vertex_count);
        vertex_count += len(vts);

    vertices = np.vstack(vertices);
    voxels = np.vstack(voxels);
    vertices, voxels = merge_identical_vertices(vertices, voxels);

    tet_mesh = form_mesh(vertices, np.array([]), voxels);
    save_mesh(tet_file, tet_mesh);

def parse_args():
    parser = argparse.ArgumentParser(description="Convert hex mesh into tet mesh");
    parser.add_argument("-s", "--symmetric", help="symmetric tet connectivity",
            action="store_true");
    parser.add_argument("hex_mesh", help="input hex mesh");
    parser.add_argument("tet_mesh", help="output tet mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    hex2tet(args.hex_mesh, args.tet_mesh, args.symmetric);

if __name__ == "__main__":
    main();

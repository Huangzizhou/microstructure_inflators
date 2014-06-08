#!/usr/bin/env python

import argparse
import numpy as np
import LinearElasticitySettings
from BoxMeshGenerator import generate_box_mesh

import PyMesh

def merge_identical_vertices(vertices, voxels):
    grid = PyMesh.HashGrid.create(1e-6);
    inverse_map = np.arange(len(vertices), dtype=int);
    unique_vertices = [];
    for i,v in enumerate(vertices):
        nearby_id = grid.get_items_near_point(v);
        if len(nearby_id) == 0:
            v_id = len(unique_vertices);
            grid.insert(v_id, v);
            inverse_map[i] = v_id;
            unique_vertices.append(v);
        else:
            inverse_map[i] = nearby_id[0];

    vertices = np.vstack(unique_vertices);
    voxels = np.vstack([inverse_map[voxel] for voxel in voxels]);
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

def hex2tet(hex_file, tet_file, num_subdiv):
    num_samples = 2**num_subdiv;
    hex_mesh = load_mesh(hex_file);
    box_meshes = [];
    hex_vertices = hex_mesh.get_vertices().reshape(
            (hex_mesh.get_num_vertices(), -1));
    hexes = hex_mesh.get_voxels().reshape((-1, hex_mesh.get_vertex_per_voxel()));

    for voxel in hexes:
        vts = hex_vertices[voxel];
        bbox_min = np.amin(vts, axis=0);
        bbox_max = np.amax(vts, axis=0);
        box_mesh = generate_box_mesh(bbox_min, bbox_max, num_samples);
        box_meshes.append(box_mesh);

    vertex_count = 0;
    vertices = [];
    voxels = [];
    for box_mesh in box_meshes:
        vertices.append(box_mesh.vertices);
        voxels.append(box_mesh.voxels + vertex_count);
        vertex_count += box_mesh.num_vertices;

    vertices = np.vstack(vertices);
    voxels = np.vstack(voxels);
    vertices, voxels = merge_identical_vertices(vertices, voxels);

    tet_mesh = form_mesh(vertices, np.array([]), voxels);
    save_mesh(tet_file, tet_mesh);

def parse_args():
    parser = argparse.ArgumentParser(description="Convert hex mesh into tet mesh");
    parser.add_argument("--subdiv", "-s", help="number of subdivision to apply",
            default=0, type=int);
    parser.add_argument("hex_mesh", help="input hex mesh");
    parser.add_argument("tet_mesh", help="output tet mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    hex2tet(args.hex_mesh, args.tet_mesh, args.subdiv);

if __name__ == "__main__":
    main();

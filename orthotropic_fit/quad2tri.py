#!/usr/bin/env python

import argparse
import numpy as np
import LinearElasticitySettings
from BoxMeshGenerator import split_quad_into_tris, split_quad_into_tris_symmetrically
from BoxMeshGenerator import remove_duplicated_vertices, remove_isolated_vertices
from BoxMeshGenerator import subdivide_quad

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

def quad2tri(quad_file, tri_file, keep_symmetry, subdiv_order):
    quad_mesh = load_mesh(quad_file);
    quad_vertices = quad_mesh.get_vertices().reshape(
            (quad_mesh.get_num_vertices(), -1));
    quads = quad_mesh.get_faces().reshape((-1, quad_mesh.get_vertex_per_face()));

    num_quads = len(quads);
    attr_names = []
    ori_attributes = [];

    for name in quad_mesh.get_attribute_names():
        attributes = quad_mesh.get_attribute(name).ravel();
        if len(attributes) % num_quads == 0:
            attr_names.append(name);
            ori_attributes.append(attributes.reshape((num_quads, -1)));

    vertex_count = 0;
    vertices = [];
    faces = [];
    quad_indices = [];
    attributes = {};

    for i, face in enumerate(quads):
        corners = quad_vertices[face];
        subcell_corners = subdivide_quad(corners, subdiv_order);
        for corners in subcell_corners:
            if keep_symmetry:
                vts, tris = split_quad_into_tris_symmetrically(corners);
            else:
                vts, tris = split_quad_into_tris(corners);
            vertices.append(vts);
            faces.append(tris + vertex_count);
            vertex_count += len(vts);
            quad_indices += [i] * len(tris);

            for name, attr in zip(attr_names, ori_attributes):
                if name not in attributes:
                    attributes[name] = [];
                attributes[name].append([attr[i]] * len(tris));

    vertices = np.vstack(vertices);
    faces = np.vstack(faces);
    vertices, faces = merge_identical_vertices(vertices, faces);
    quad_indices = np.array(quad_indices);

    tri_mesh = form_mesh(vertices, faces, np.array([]));
    tri_mesh.add_attribute("cell_index");
    tri_mesh.set_attribute("cell_index", quad_indices);

    for attr_name in attributes:
        attr = np.vstack(attributes[attr_name]).ravel(order="C");
        tri_mesh.add_attribute(attr_name);
        tri_mesh.set_attribute(attr_name, attr);

    save_mesh(tri_file, tri_mesh, *tri_mesh.get_attribute_names());

def parse_args():
    parser = argparse.ArgumentParser(description="Convert quad mesh into triangle mesh");
    parser.add_argument("-s", "--symmetric", help="symmetric triangle connectivity",
            action="store_true");
    parser.add_argument("--subdiv", help="subdivision order", type=int, default=0);
    parser.add_argument("quad_mesh", help="input quad mesh");
    parser.add_argument("tri_mesh", help="output tri mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    quad2tri(args.quad_mesh, args.tri_mesh, args.symmetric, args.subdiv);

if __name__ == "__main__":
    main();

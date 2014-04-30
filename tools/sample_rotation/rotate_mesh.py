#!/usr/bin/env python

import argparse
from math import sin, cos, radians
import numpy as np
from numpy.linalg import norm

import PyMeshSetting
import PyMesh
from quaternion import Quaternion

def load_mesh(mesh_file):
    factory = PyMesh.MeshFactory();
    factory.load_file(mesh_file);
    mesh = factory.create();
    return mesh;

def form_mesh(vertices, faces, voxels=np.array([])):
    factory = PyMesh.MeshFactory();
    factory.load_data(
            vertices.ravel(order="C"),
            faces.ravel(order="C"),
            voxels.ravel(order="C"),
            3, 3, 4);
    mesh = factory.create();
    return mesh;

def save_mesh(mesh, output_file):
    writer = PyMesh.MeshWriter.create_writer(output_file);
    writer.write_mesh(mesh);

def form_rotation_matrix(theta, phi):
    theta = radians(theta);
    phi = radians(phi);
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);
    v = np.array([x, y, z]);
    v = v / norm(v);
    z_axis = np.array([0.0, 0.0, 1.0]);
    quaternion = Quaternion.fromData(z_axis, v);
    return quaternion.to_matrix();

def rotate(mesh, theta, phi):
    dim = mesh.get_dim();
    vertices = mesh.get_vertices().reshape((-1, dim), order="C");
    rot_mat = form_rotation_matrix(theta, phi);
    vertices = np.dot(rot_mat, vertices.T).T;
    rot_mesh = form_mesh(vertices, mesh.get_faces());
    return rot_mesh;

def parse_args():
    parser = argparse.ArgumentParser(description="Generate rotated mesh");
    parser.add_argument("input", help="input model");
    parser.add_argument("output", help="output model");
    parser.add_argument("--theta", default=0.0, type=float,
            help="spherical coordinate -- z angle");
    parser.add_argument("--phi", default=0.0, type=float,
            help="spherical coordinate -- xy angle");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.input);
    mesh = rotate(mesh, args.theta, args.phi);
    save_mesh(mesh, args.output);

if __name__ == "__main__":
    main();

#!/usr/bin/env python

import sys
import os
import numpy as np

sys.path.append("../lib");
sys.path.append("../swig");

import PyWireInflator2D
import PyMeshSetting
import PyMesh

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

def main():
    inflator = PyWireInflator2D.WireInflatorFacade();
    num_params = inflator.get_num_parameters();
    param = np.zeros(num_params);

    # Thickness
    param[:5] = 0.05;
    param[0] = 0.10;
    #param[4] = 0.05;
    #param[3] = 0.05;
    # Offset
    param[5:] = 0.0;
    param[7] = 0.1;

    inflator.set_parameter(param);
    inflator.generate_periodic_pattern();

    vertices = inflator.get_vertices().reshape((-1, 2));
    triangles = inflator.get_triangles().reshape((-1, 3));
    velocity = inflator.get_boundary_velocity();

    attribute_names = [];
    mesh = form_mesh(vertices, triangles);
    mesh.add_attribute("vertex_normal");
    normals = mesh.get_attribute("vertex_normal").reshape((-1, 2), order="C");
    for i,v in enumerate(velocity.T):
        name = "velocity_{}".format(i);
        mesh.add_attribute(name);
        grad = normals * v[:, np.newaxis];
        mesh.set_attribute(name, grad.ravel());
        attribute_names.append(name);

    save_mesh("tmp.msh", mesh, *attribute_names);

if __name__ == "__main__":
    main();

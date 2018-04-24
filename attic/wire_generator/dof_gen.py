#!/usr/bin/env python

import argparse
import numpy as np
import os.path

from tile import load_wire, load_mesh
from parameter.PyParameters import PyParameters

import PyMesh
import PyWires

def load_parameters(wire_network, dof_type, thickness_type):
    parameters = PyParameters(wire_network, 0.5);
    if thickness_type == "vertex":
        thickness_type = PyWires.VERTEX;
    else:
        thickness_type = PyWires.EDGE;

    if dof_type == "orthotropic":
        parameters.load_default_orthotropic_parameters(thickness_type);
    else:
        parameters.load_default_isotropic_parameters(thickness_type);

    return parameters;

def save_mesh(mesh, mesh_file):
    basename, ext = os.path.splitext(mesh_file);
    writer = PyMesh.MeshWriter.create_writer(mesh_file);
    if ext in (".msh", ".ply"):
        attribute_names = mesh.get_attribute_names();
        for attr_name in attribute_names:
            writer.with_attribute(attr_name);

    writer.in_ascii();
    writer.write_mesh(mesh);

def parse_args():
    parser = argparse.ArgumentParser(
            description="Generate dof attribute for a given mesh");
    parser.add_argument("--dof-type", choices=["isotropic", "orthotropic"],
            help="type of dof", default="orthotropic");
    parser.add_argument("--thickness-type", choices=["vertex", "edge"],
            help="per-vertex or per-edge thickness", default="vertex");
    parser.add_argument("wire_network", help="target wire network");
    parser.add_argument("guide_mesh", help="guide mesh");
    parser.add_argument("out_mesh", help="output mesh");
    return parser.parse_args();

def main():
    args = parse_args();

    wire_network = load_wire(args.wire_network);
    guide_mesh = load_mesh(args.guide_mesh);
    parameters = load_parameters(wire_network, args.dof_type,
            args.thickness_type);

    if guide_mesh.get_dim() == 2:
        num_elems = guide_mesh.get_num_faces();
    else:
        num_elems = guide_mesh.get_num_voxels();

    dofs = np.array([parameters.dofs] * num_elems);
    for i in range(parameters.num_dofs):
        dof_name = "dof_{}".format(i);
        guide_mesh.add_attribute(dof_name);
        guide_mesh.set_attribute(dof_name, dofs[:, i].ravel());

    save_mesh(guide_mesh, args.out_mesh);

if __name__ == "__main__":
    main();

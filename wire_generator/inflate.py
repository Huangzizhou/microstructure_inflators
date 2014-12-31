#!/usr/bin/env python
import argparse
import numpy as np
import os.path

from core.WireNetwork import WireNetwork
from inflator.WireInflator import WireInflator
from inflator.PyWiresInflator import PyWiresInflator
from parameter.PyParameters import PyParameters

import PyMesh

def save_mesh(mesh_file, mesh):
    basename, ext = os.path.splitext(mesh_file);
    writer = PyMesh.MeshWriter.create_writer(mesh_file);
    if ext in (".msh", ".ply"):
        attr_names = mesh.get_attribute_names();
        for name in attr_names:
            writer.with_attribute(name);
    writer.write_mesh(mesh);

def inflate(wire_file, thickness, mesh_file, with_symmetry_orbits, subdiv_order,
        subdiv_method):
    wires = WireNetwork();
    wires.load_from_file(wire_file);
    wires.add_attribute("thickness", np.ones(wires.num_vertices) * thickness);
    parameters = PyParameters(wires, thickness);
    inflator = PyWiresInflator(wires, parameters);
    inflator.inflate(True, subdiv_order, subdiv_method);
    mesh = inflator.mesh;
    save_mesh(mesh_file, mesh);

def parse_args():
    parser = argparse.ArgumentParser(description="Inflate a wire network");
    parser.add_argument("--thickness", "-t", help="target wire thickness",
            required=True, type=float);
    parser.add_argument("--output", "-o", help="output mesh file", required=True);
    parser.add_argument("--subdiv", "-s", type=int, default=0,
            help="number of subdivisions");
    parser.add_argument("--subdiv-method", help="subdivision mehtod",
            choices=["simple", "loop"], default="simple");
    parser.add_argument("--with-symmetry-orbits",
            help="output also symmetry orbits", action="store_true");
    parser.add_argument("wire_file", help="input .wire file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    inflate(args.wire_file, args.thickness, args.output,
            args.with_symmetry_orbits, args.subdiv, args.subdiv_method);

if __name__ == "__main__":
    main();

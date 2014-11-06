#!/usr/bin/env python
import argparse
import numpy as np

from core.WireNetwork import WireNetwork
from inflator.WireInflator import WireInflator

def inflate(wire_file, thickness, mesh_file, with_symmetry_orbits, subdiv_order):
    wires = WireNetwork();
    wires.load_from_file(wire_file);
    wires.attributes["vertex_thickness"] = \
            np.ones(wires.num_vertices) * thickness;
    if with_symmetry_orbits:
        wires.attributes.add("symmetry_vertex_orbit");
        wires.attributes.add("symmetry_edge_orbit");
    inflator = WireInflator(wires);
    inflator.inflate(thickness, subdiv_order, "loop");
    inflator.save(mesh_file);

def parse_args():
    parser = argparse.ArgumentParser(description="Inflate a wire network");
    parser.add_argument("--thickness", "-t", help="target wire thickness",
            required=True, type=float);
    parser.add_argument("--output", "-o", help="output mesh file", required=True);
    parser.add_argument("--subdiv", "-s", type=int, default=0,
            help="number of subdivisions");
    parser.add_argument("--with-symmetry-orbits",
            help="output also symmetry orbits", action="store_true");
    parser.add_argument("wire_file", help="input .wire file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    inflate(args.wire_file, args.thickness, args.output,
            args.with_symmetry_orbits, args.subdiv);

if __name__ == "__main__":
    main();

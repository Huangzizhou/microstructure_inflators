#!/usr/bin/env python
import argparse

from WireNetwork import WireNetwork
from WireInflator import WireInflator

def inflate(wire_file, thickness, mesh_file):
    wires = WireNetwork();
    wire.load_from_file(wire_file);
    inflator = WireInflator(wires);
    inflator.inflate(thickness);
    inflator.save(mesh_file);

def parse_args():
    parser = argparse.ArgumentParser(description="Inflate a wire network");
    parser.add_argument("--thickness", "-t", help="target wire thickness",
            required=True, type=float);
    parser.add_argument("--output", "-o", help="output mesh file", required=True);
    parser.add_argument("wire_file", help="input .wire file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    inflate(args.wire_file, args.thickness, args.output);

if __name__ == "__main__":
    main();

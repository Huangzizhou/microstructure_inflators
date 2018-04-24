#!/usr/bin/env python3.6

# System libs
import os
import argparse
import subprocess

# Local libs
import paths

def convert(input_mesh, output_mesh):
    converter = paths.find("mesh_convert", paths.MESHFEM_DIR)
    if converter is None:
        raise FileNotFoundError("mesh_convert")

    args = [converter, input_mesh, output_mesh]
    subprocess.check_call(args)

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="input mesh")
    parser.add_argument("output", help="output mesh")
    return parser.parse_args()


def main():
    args = parse_args()
    convert(args.input, args.output)

if __name__ == "__main__":
    main()
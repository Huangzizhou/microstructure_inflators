#!/usr/bin/env python

import argparse
import os
from subprocess import check_call

PROJECT_DIR = os.environ["MICROSTRUCTURES_PATH"];

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate material hex mesh");
    parser.add_argument("--output", "-o", help="output hex mesh");
    parser.add_argument("--young", help="young attribute name", default="Final E");
    parser.add_argument("--poisson", help="poisson attribute name",
            default="Final nu");
    parser.add_argument("material_mesh",
            help="tet mesh containing target material properties");
    parser.add_argument("hex_mesh",
            help="hex mesh that defines cells");
    return parser.parse_args();

def main():
    args = parse_args();
    command = "cp {} {}".format(args.hex_mesh, args.output);
    check_call(command.split());

    exe_name = os.path.join(PROJECT_DIR, "tools/msh_tools/rename_element_attribute.py");
    command = [exe_name, args.material_mesh, args.young, "young"];
    check_call(command);
    command = [exe_name, args.material_mesh, args.poisson, "poisson"];
    check_call(command);

    exe_name = os.path.join(PROJECT_DIR,
            "tools/msh_tools/move_element_attribute.py");
    command = "{} --attribute young --attribute poisson {} {}".format(
            exe_name, args.material_mesh, args.output);
    check_call(command.split());

if __name__ == "__main__":
    main();

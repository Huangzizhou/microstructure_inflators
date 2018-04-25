#!/usr/bin/env python

import argparse
import numpy as np
import parser
import re

import PyMeshSetting
import PyMesh

from add_element_attribute import load_mesh, save_mesh

def parse_args():
    parser = argparse.ArgumentParser(
            description="Rename attribute");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("old_attr_name", help="old attribute name");
    parser.add_argument("new_attr_name", help="new attribute name");
    return parser.parse_args();

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);

    assert(mesh.has_attribute(args.old_attr_name));
    attr_value = mesh.get_attribute(args.old_attr_name);
    mesh.add_attribute(args.new_attr_name);
    mesh.set_attribute(args.new_attr_name, attr_value);
    #mesh.remove_attribute(args.old_attr_name);

    save_mesh(args.input_mesh, mesh, *mesh.get_attribute_names());

if __name__ == "__main__":
    main();

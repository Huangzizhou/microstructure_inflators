#!/usr/bin/env python

import argparse
import numpy as np

import PyMeshSetting
import PyMesh

from add_element_attribute import load_mesh, save_mesh

def is_hex(mesh):
    return mesh.get_vertex_per_voxel() == 8;

def is_tet(mesh):
    return mesh.get_vertex_per_voxel() == 4;

def is_quad(mesh):
    return mesh.get_vertex_per_face() == 4;

def is_tri(mesh):
    return mesh.get_vertex_per_face() == 3;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Move attribute from input mesh to output mesh");
    parser.add_argument("--attribute", "-a", help="attribute to move",
            action="append");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    return parser.parse_args();

def main():
    args = parse_args();
    in_mesh = load_mesh(args.input_mesh);
    out_mesh = load_mesh(args.output_mesh);

    for attr_name in args.attribute:
        assert(in_mesh.has_attribute(attr_name));
        attr_value = in_mesh.get_attribute(attr_name);

        if is_tet(in_mesh) and is_tet(out_mesh):
            out_mesh.add_attribute(attr_name);
            out_mesh.set_attribute(attr_name, attr_value);
        elif is_tri(in_mesh) and is_tri(out_mesh):
            out_mesh.add_attribute(attr_name);
            out_mesh.set_attribute(attr_name, attr_value);
        elif is_tet(in_mesh) and is_hex(out_mesh):
            assert(in_mesh.has_attribute("cell_index"));
            index_map = in_mesh.get_attribute("cell_index").ravel().astype(int);

            attr_value = attr_value.reshape((in_mesh.get_num_voxels(), -1),
                    order="C");
            hex_attr_value = np.zeros((out_mesh.get_num_voxels(),
                attr_value.shape[1]));
            hex_attr_value[index_map] = attr_value;

            out_mesh.add_attribute(attr_name);
            out_mesh.set_attribute(attr_name, hex_attr_value);
        elif is_tri(in_mesh) and is_quad(out_mesh):
            assert(in_mesh.has_attribute("cell_index"));
            index_map = in_mesh.get_attribute("cell_index").ravel().astype(int);

            attr_value = attr_value.reshape((in_mesh.get_num_faces(), -1),
                    order="C");
            hex_attr_value = np.zeros((out_mesh.get_num_faces(),
                attr_value.shape[1]));
            hex_attr_value[index_map] = attr_value;

            out_mesh.add_attribute(attr_name);
            out_mesh.set_attribute(attr_name, hex_attr_value);
        else:
            print(in_mesh.get_vertex_per_face());
            print(out_mesh.get_vertex_per_face());
            assert(False);

    save_mesh(args.output_mesh, out_mesh, *out_mesh.get_attribute_names());

if __name__ == "__main__":
    main();

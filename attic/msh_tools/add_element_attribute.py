#!/usr/bin/env python

import argparse
import numpy as np
import parser
import re

import PyMeshSetting
import PyMesh

def load_mesh(mesh_file):
    factory = PyMesh.MeshFactory();
    factory.load_file(mesh_file);
    factory.drop_zero_dim();
    return factory.create();

def save_mesh(mesh_file, mesh, *attributes):
    writer = PyMesh.MeshWriter.create_writer(mesh_file);
    for attr in attributes:
        assert(mesh.has_attribute(attr));
        writer.with_attribute(attr);
    writer.write_mesh(mesh);

class Formula:
    def __init__(self, mesh):
        self.mesh = mesh;
        self.field_parser = re.compile("\{\w+\}");
        if self.mesh.get_dim() == 2:
            self.index = np.arange(self.mesh.get_num_faces(), dtype=int);
        elif self.mesh.get_dim() == 3:
            self.index = np.arange(self.mesh.get_num_voxels(), dtype=int);

    def eval_formula(self, formula):
        self.formula = formula;
        field_names = self.field_parser.findall(formula);
        for name in field_names:
            field_name = name[1:-1];
            if hasattr(self, field_name):
                escaped_name = field_name;
            else:
                escaped_name = self.load_field(field_name);
            self.formula = self.formula.replace(name,
                    "self.{}".format(escaped_name));

        st = parser.expr(self.formula);
        code = st.compile("file.py");
        return eval(code);

    def load_field(self, name):
        assert(self.mesh.has_attribute(name));
        value = self.mesh.get_attribute(name);
        escaped_name = name.replace(" ", "_")
        setattr(self, name, value);
        return escaped_name

def parse_args():
    parser = argparse.ArgumentParser(description="Added element attribute to mesh");
    parser.add_argument("--name", help="attribute name");
    parser.add_argument("--formula", help="formula to compute new attribute");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    return parser.parse_args();

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);

    formula = Formula(mesh);
    values = formula.eval_formula(args.formula);
    if isinstance(values, list):
        values = np.array(values);
    elif isinstance(values, (int, float)):
        values = np.ones(mesh.get_num_voxels()) * values;
    mesh.add_attribute(args.name);
    mesh.set_attribute(args.name, values);

    save_mesh(args.output_mesh, mesh, *mesh.get_attribute_names());

if __name__ == "__main__":
    main();

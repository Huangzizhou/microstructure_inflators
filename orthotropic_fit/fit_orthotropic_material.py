#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from OrthotropicFitter import OrthotropicFitter
from timethis import timethis

def fit_orthotropic(mesh):
    fitter = OrthotropicFitter(mesh);
    fitter.fit();
    pressures = fitter.pressures;
    displacements = fitter.displacements;
    stress_traces = fitter.stress_traces;
    young = fitter.youngs_modulus;
    poisson = fitter.poisson_ratio;
    shear = fitter.shear_modulus;
    error = fitter.residual_error;
    condition_num = fitter.condition_num;
    return pressures, displacements, stress_traces,\
            young, poisson, shear, error, condition_num;

def add_attribute_to_mesh(mesh, attr_name, attr_value):
    if mesh.has_attribute(attr_name):
        raise RuntimeError("Attribute {} already exists.".format(attr_name));
    mesh.add_attribute(attr_name);
    mesh.set_attribute(attr_name, attr_value);

def save_mesh_fields(mesh_file, mesh, pressures, displacements, stress_traces):
    num_fields = len(displacements);
    attributes_to_save = [];
    for i in range(num_fields):
        pressure_attr_name = "pressure_{}".format(i);
        disp_attr_name = "displacement_{}".format(i);
        strs_attr_name = "stress_trace_{}".format(i);
        add_attribute_to_mesh(mesh, pressure_attr_name, pressures[i]);
        add_attribute_to_mesh(mesh, disp_attr_name, displacements[i]);
        add_attribute_to_mesh(mesh, strs_attr_name, stress_traces[i]);
        attributes_to_save.append(pressure_attr_name);
        attributes_to_save.append(disp_attr_name);
        attributes_to_save.append(strs_attr_name);
    save_mesh(mesh_file, mesh, *attributes_to_save);

def save_parameters(mesh_file, young, poisson, shear, error, condition_num):
    basename, ext = os.path.splitext(mesh_file);
    parameter_file = basename + ".param";
    parameter_dict = {
            "youngs_modulus": young,
            "poisson_ratio": poisson,
            "shear_modulus": shear,
            "residual_error": error,
            "condition_num": condition_num }
    with open(parameter_file, 'w') as fout:
        json.dump(parameter_dict, fout, indent=4);

def parse_args():
    parser = argparse.ArgumentParser(
            description="Compute the best fit orthotropic material parameters");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);
    pressures, displacements, stress_traces,\
            young, poisson, shear, error, condition_num = fit_orthotropic(mesh)
    save_mesh_fields(args.output_mesh, mesh, pressures, displacements, stress_traces);
    save_parameters(args.output_mesh,
            young, poisson, shear, error, condition_num);

if __name__ == "__main__":
    main();
    timethis.summarize();

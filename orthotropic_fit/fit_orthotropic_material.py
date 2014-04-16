#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

from HomogenizationValidator import HomogenizationValidator
import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from OrthotropicFitter import OrthotropicFitter
from timethis import timethis

import PyAssembler

@timethis
def fit_orthotropic_material_and_validate(input_file, output_file, material_file):
    mesh = load_mesh(input_file);
    fitter = fit(mesh, material_file);
    validator = validate(fitter);

    basename,ext = os.path.splitext(output_file);
    validation_file = basename + "_homogenized" + ext;
    param_file = basename + "_param.json";
    coarse_mesh_file = basename + "_coarse" + ext;

    save_mesh_fields(output_file, fitter.mesh,
            fitter.pressures,
            fitter.displacements,
            fitter.stress_traces);

    #save_mesh_fields(coarse_mesh_file, fitter.coarse_mesh,
    #        displacements = fitter.coarse_displacements);

    save_mesh_fields(validation_file, validator.mesh,
            validator.pressures,
            validator.displacements,
            validator.stress_traces);

    save_parameters(param_file,
            fitter.youngs_modulus.tolist(),
            fitter.poisson_ratio.tolist(),
            fitter.shear_modulus.tolist(),
            fitter.residual_error,
            fitter.condition_num);

@timethis
def validate(fitter):
    mat = PyAssembler.Material.create_orthotropic(1.0,
            fitter.youngs_modulus,
            fitter.poisson_ratio,
            fitter.shear_modulus);
    validator = HomogenizationValidator(fitter.mesh, mat);
    validator.simulate(fitter.bc_configs);
    return validator;

@timethis
def fit(mesh, material_file):
    fitter = OrthotropicFitter(mesh, material_file);
    fitter.fit();
    return fitter;

@timethis
def add_attribute_to_mesh(mesh, attr_name, attr_value):
    if mesh.has_attribute(attr_name):
        raise RuntimeError("Attribute {} already exists.".format(attr_name));
    mesh.add_attribute(attr_name);
    mesh.set_attribute(attr_name, attr_value);

@timethis
def save_mesh_fields(mesh_file, mesh, pressures=None, displacements=None,
        stress_traces=None):
    num_fields = len(displacements);
    attributes_to_save = [];
    for i in range(num_fields):
        if pressures is not None:
            pressure_attr_name = "pressure_{}".format(i);
            add_attribute_to_mesh(mesh, pressure_attr_name, pressures[i]);
            attributes_to_save.append(pressure_attr_name);

        if displacements is not None:
            disp_attr_name = "displacement_{}".format(i);
            add_attribute_to_mesh(mesh, disp_attr_name, displacements[i]);
            attributes_to_save.append(disp_attr_name);

        if stress_traces is not None:
            strs_attr_name = "stress_trace_{}".format(i);
            add_attribute_to_mesh(mesh, strs_attr_name, stress_traces[i]);
            attributes_to_save.append(strs_attr_name);
    save_mesh(mesh_file, mesh, *attributes_to_save);

@timethis
def save_parameters(parameter_file, young, poisson, shear, error, condition_num):
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
    parser.add_argument("--material", help="material file", default=None);
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    fit_orthotropic_material_and_validate(args.input_mesh, args.output_mesh,
            args.material);

if __name__ == "__main__":
    main();
    timethis.summarize();

#!/usr/bin/env python

import argparse
import numpy as np
import os.path

from HomogenizationValidator import HomogenizationValidator
import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from MaterialFitterFactory import MaterialFitterFactory
from ResultOutputUtils import save_mesh_fields, save_parameters
from timethis import timethis

import PyAssembler

@timethis
def fit_material_and_validate(input_file, output_file, material_model, material_file):
    mesh = load_mesh(input_file);
    fitter = fit(mesh, material_model, material_file);

    basename,ext = os.path.splitext(output_file);
    validation_file = basename + "_homogenized" + ext;
    param_file = basename + "_param.json";
    coarse_mesh_file = basename + "_coarse" + ext;

    save_mesh_fields(output_file, fitter.mesh,
            fitter.pressures,
            fitter.displacements,
            fitter.stress_traces);

    save_mesh_fields(coarse_mesh_file, fitter.coarse_mesh,
            displacements = fitter.coarse_displacements);

    validator = validate(fitter, material_model);
    save_mesh_fields(validation_file, validator.mesh,
            validator.pressures,
            validator.displacements,
            validator.stress_traces);

    save_parameters(param_file, fitter);

@timethis
def validate(fitter, material_model):
    if material_model == "orthotropic":
        mat = PyAssembler.Material.create_orthotropic(1.0,
                fitter.youngs_modulus,
                fitter.poisson_ratio,
                fitter.shear_modulus);
    elif material_model == "symmetric":
        mat = PyAssembler.Material.create_symmetric(1.0,
                fitter.elasticity_tensor);
    validator = HomogenizationValidator(fitter.mesh, mat);
    validator.simulate(fitter.bc_configs);
    return validator;

@timethis
def fit(mesh, material_model, material_file):
    factory = MaterialFitterFactory(mesh, material_file);
    fitter = factory.create(material_model);
    fitter.fit();
    return fitter;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Compute the best fit material parameters");
    parser.add_argument("--material", help="base material file", default=None);
    parser.add_argument("--material-model", help="target material model",
            choices=["orthotropic", "symmetric"], default="orthotropic");
    parser.add_argument("input_mesh", help="input mesh");
    parser.add_argument("output_mesh", help="output mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    fit_material_and_validate(args.input_mesh, args.output_mesh,
            args.material_model, args.material);

if __name__ == "__main__":
    main();
    timethis.summarize();

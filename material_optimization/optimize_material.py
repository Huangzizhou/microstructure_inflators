#!/usr/bin/env python

import argparse
from scipy.optimize import minimize
import numpy as np
from numpy.linalg import norm

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from BoundaryConditionExtractor import BoundaryConditionExtractor
from timethis import timethis

from IsotropicMatOptSetting import IsotropicMatOptSetting

def optimize(setting, max_iterations):
    init_parameters = setting.parameters;
    result = minimize(setting.evaluate, init_parameters, method="BFGS", jac=True,
            callback=setting.log_iteration,
            options={
                "maxiter": max_iterations,
                "disp": True
                });
    return result.x;

def create_setting(mesh, bd_file):
    bc_extractor = BoundaryConditionExtractor(mesh);
    bc_extractor.extract_from_file(bd_file);

    setting = IsotropicMatOptSetting(mesh,
            bc_extractor.neumann_bc, bc_extractor.dirichlet_bc);
    return setting;

def save_result(optimizer, out_mesh_name, *attribute_names):
    for name in attribute_names:
        if not optimizer.mesh.has_attribute(name):
            optimizer.mesh.add_attribute(name);
            optimizer.mesh.set_attribute(name, getattr(optimizer, name));

    save_mesh(out_mesh_name, optimizer.mesh,
            *attribute_names);

def dump_iteration_log(setting, log_file):
    import csv
    with open(log_file, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(["iteration", "objective", "gradient"]);

        for i,idx in enumerate(setting.iteration_indices):
            writer.writerow([i,
                setting.objective_history[idx],
                norm(setting.gradient_history[idx]) ]);

def parse_args():
    parser = argparse.ArgumentParser(description="Optimize material properties");
    parser.add_argument("-n", "--num-iterations", default=10, type=int);
    parser.add_argument("-b", "--boundary-condition",
            help="Boundary condition specification", required=True);
    parser.add_argument("--log-file", help="log file name",
            default="iteration_history.csv");
    parser.add_argument("input_mesh", help="input mesh file");
    parser.add_argument("output_mesh", help="output mesh file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);
    setting = create_setting(mesh, args.boundary_condition);
    opt_solution = optimize(setting, args.num_iterations);

    attribute_names = [];
    num_elements = mesh.num_elements;
    for i,idx in enumerate(setting.iteration_indices):
        young = setting.parameter_history[idx][:num_elements];
        young_name = "young_{}".format(i);
        mesh.add_attribute(young_name);
        mesh.set_attribute(young_name, young);
        attribute_names.append(young_name);

        poisson = setting.parameter_history[idx][num_elements:];
        poisson_name = "poisson_{}".format(i);
        mesh.add_attribute(poisson_name);
        mesh.set_attribute(poisson_name, poisson);
        attribute_names.append(poisson_name);

    save_result(setting, args.output_mesh,
            "young", "poisson",
            "displacement", "source_term",
            "lagrange_multiplier", "lagrange_source_term",
            "grad_young", "grad_poisson",
            "target_displacement",
            *attribute_names);

    dump_iteration_log(setting, args.log_file);


if __name__ == "__main__":
    main();
    timethis.summarize();

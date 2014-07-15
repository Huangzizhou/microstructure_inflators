#!/usr/bin/env python

import argparse
from scipy.optimize import minimize
import numpy as np
from numpy.linalg import norm
import os.path

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from BoundaryConditionExtractor import BoundaryConditionExtractor
from timethis import timethis

from IsotropicMatOptSetting import IsotropicMatOptSetting

def optimize(setting, max_iterations):
    GTOL = 1e-2;
    init_parameters = setting.parameters;
    result = minimize(setting.evaluate, init_parameters,
            method="L-BFGS-B",
            jac=True,
            bounds = setting.bounds,
            callback=setting.log_iteration,
            options={
                "maxiter": max_iterations,
                "gtol": min(1e-6, GTOL / setting.mesh.num_elements),
                #"gtol": GTOL / setting.mesh.num_elements,
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

def dump_iteration_log(setting, output_file):
    basename, ext = os.path.splitext(output_file);
    log_file = basename + ".csv";
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

    def add_output_field(name, value):
        mesh.add_attribute(name);
        mesh.set_attribute(name, value);
        attribute_names.append(name);

    for i,idx in enumerate(setting.iteration_indices):
        young = setting.parameter_history[idx][:num_elements];
        young_name = "young_{}".format(i);
        add_output_field(young_name, young);

        poisson = setting.parameter_history[idx][num_elements:];
        poisson_name = "poisson_{}".format(i);
        add_output_field(poisson_name, poisson);

        displacement = setting.displacement_history[idx];
        displacement_name = "disp_{}".format(i)
        add_output_field(displacement_name, displacement);

        grad_young = setting.grad_young_history[idx];
        grad_young_name = "grad_young_{}".format(i);
        add_output_field(grad_young_name, grad_young);

        lagrange_multiplier = setting.lagrange_history[idx];
        lagrange_name = "lagrange_{}".format(i);
        add_output_field(lagrange_name, lagrange_multiplier);

        displacement_strain = setting.displacement_strain_history[idx];
        displacement_strain_name = "u_strain_{}".format(i);
        add_output_field(displacement_strain_name, displacement_strain);

        lagrange_strain = setting.lagrange_strain_history[idx];
        lagrange_strain_name = "l_strain_{}".format(i);
        add_output_field(lagrange_strain_name, lagrange_strain);

    save_result(setting, args.output_mesh,
            "young", "poisson",
            "displacement", "source_term",
            "lagrange_multiplier", "lagrange_source_term",
            "grad_young", "grad_poisson",
            "target_displacement",
            *attribute_names);

    dump_iteration_log(setting, args.output_mesh);


if __name__ == "__main__":
    main();
    timethis.summarize();

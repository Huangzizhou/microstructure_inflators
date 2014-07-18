#!/usr/bin/env python

import argparse
import csv
import os.path
from MaterialOptimizer import MaterialOptimizer

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from BoundaryConditionExtractor import BoundaryConditionExtractor

def add_attribute(mesh, name, value, attr_names):
    if not mesh.has_attribute(name):
        mesh.add_attribute(name);
    mesh.set_attribute(name, value);
    attr_names.append(name);

def save_results(optimizer, output_file):
    attr_names = [];
    mesh = optimizer.mesh;
    add_attribute(mesh, "primal displacement", optimizer.primal.displacement,
            attr_names);
    add_attribute(mesh, "dual displacement", optimizer.dual.displacement,
            attr_names);
    add_attribute(mesh, "source term", optimizer.dual.source_term, attr_names);

    for i in range(optimizer.iterations):
        add_attribute(mesh, "young_{:03}".format(i), optimizer.young[i],
                attr_names);
        add_attribute(mesh, "poisson_{:03}".format(i), optimizer.poisson[i],
                attr_names);
        add_attribute(mesh, "primal_u_{:03}".format(i),
                optimizer.primal_displacement[i], attr_names);
        add_attribute(mesh, "dual_u_{:03}".format(i),
                optimizer.dual_displacement[i], attr_names);

    save_mesh(output_file, mesh,
            optimizer.material.young_attr_name,
            optimizer.material.poisson_attr_name,
            *attr_names);

def save_progress(optimizer, output_file):
    basename, ext = os.path.splitext(output_file);
    csv_file = basename + ".csv";
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(["iteration", "objective"]);
        for i,obj in enumerate(optimizer.objective_history):
            writer.writerow([i, obj]);

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
    bc_extractor = BoundaryConditionExtractor(mesh);
    bc_extractor.extract_from_file(args.boundary_condition);

    optimizer = MaterialOptimizer(mesh,
            bc_extractor.dirichlet_bc,
            bc_extractor.neumann_bc);

    optimizer.optimize(args.num_iterations);
    save_results(optimizer, args.output_mesh);
    save_progress(optimizer, args.output_mesh);

if __name__ == "__main__":
    main();

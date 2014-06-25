#!/usr/bin/env python

import argparse
import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from BoundaryConditionExtractor import BoundaryConditionExtractor
from timethis import timethis

from IsotropicMaterialOptimizer import IsotropicMaterialOptimizer

def init_boundary_condition(optimizer, bd_file):
    bc_extractor = BoundaryConditionExtractor(optimizer.mesh);
    bc_extractor.extract_from_file(bd_file);

    optimizer.add_neumann_bc(bc_extractor.neumann_bc);
    optimizer.add_dirichlet_bc(bc_extractor.dirichlet_bc);

def save_result(optimizer, out_mesh_name, *attribute_names):
    for name in attribute_names:
        if not optimizer.mesh.has_attribute(name):
            optimizer.mesh.add_attribute(name);
            optimizer.mesh.set_attribute(name, getattr(optimizer, name));

    save_mesh(out_mesh_name, optimizer.mesh,
            *attribute_names);

def parse_args():
    parser = argparse.ArgumentParser(description="Optimize material properties");
    parser.add_argument("-b", "--boundary-condition",
            help="Boundary condition specification", required=True);
    parser.add_argument("input_mesh", help="input mesh file");
    parser.add_argument("output_mesh", help="output mesh file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.input_mesh);
    optimizer = IsotropicMaterialOptimizer(mesh);
    init_boundary_condition(optimizer, args.boundary_condition);
    optimizer.optimize(10);
    save_result(optimizer, args.output_mesh,
            "young", "poisson",
            "displacement", "source_term",
            "lagrange_multiplier", "lagrange_source_term",
            "grad_young", "grad_poisson",
            "target_displacement"
            );

if __name__ == "__main__":
    main();
    timethis.summarize();

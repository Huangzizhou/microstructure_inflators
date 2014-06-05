#!/usr/bin/env python

import sys
import argparse
import os
import os.path
import numpy as np
from numpy.linalg import norm

PYMESH_PATH = os.environ.get("PYMESH_PATH");
if PYMESH_PATH is None:
    raise ImportError("Please set PYMESH_PATH to the correct lib path.")
sys.path.append(os.path.join(PYMESH_PATH, "lib"));
sys.path.append(os.path.join(PYMESH_PATH, "swig"));
LINEAR_ELASTICITY_PATH = os.environ.get("LINEAR_ELASTICITY_PATH");
sys.path.append(LINEAR_ELASTICITY_PATH);
sys.path.append(os.path.join(LINEAR_ELASTICITY_PATH, "PyUtils"));

from BoundaryConditionExtractor import BoundaryConditionExtractor
from ElasticityUtils import displacement_to_stress, stress_norm,\
        displacement_norm, force_to_pressure
import ElasticModel2
from Material import Material
from MeshUtils import elem2vertex
from timethis import timethis
from mesh_io import load_mesh
import PyMesh
import PyAssembler
from PeriodicHomogenization import PeriodicHomogenization

def init(mesh_file, material_file):
    mesh = load_mesh(mesh_file);
    material = Material(mesh.dim, material_file);
    assembler = PyAssembler.FEAssembler.create(mesh.raw_mesh,
            material.material);
    assembler = ElasticModel2.PyAssembler(mesh, assembler, material);
    return mesh, assembler;

@timethis
def homogenize(mesh_file, material_file):
    mesh, assembler = init(mesh_file, material_file);

    homogenizer = PeriodicHomogenization(mesh, assembler);
    E = homogenizer.periodicHomogenize();
    print E
    print "Moduli: "
    print 1.0 / np.diag(np.linalg.inv(E))

def parse_args():
    parser = argparse.ArgumentParser(description="Periodic Homogenization");
    parser.add_argument("--timing", action="store_true",
            help="Display timing information");
    parser.add_argument("--material", help="Material file", default=None);
    parser.add_argument("in_mesh", help="input mesh file (.msh)");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    homogenize(args.in_mesh, args.material)

    if args.timing:
        timethis.summarize();

if __name__ == "__main__":
    main();


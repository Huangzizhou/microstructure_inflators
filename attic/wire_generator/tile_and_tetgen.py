#!/usr/bin/env python

import argparse

import core.PyMeshSetting
import PyMesh
from tile import tile, parse_config_file
from utils.mesh_io import save_mesh_raw
from utils.Tetgen import tetgen
from utils.timethis import timethis

def generate_tetgen_flags(config):
    flags = "qpQY";
    if "subdiv" in config:
        order = config["subdiv"];
        flags += "a{}".format(0.125 / (8**order));
    return flags;

@timethis
def run_tetgen(config, mesh):
    flags = generate_tetgen_flags(config);

    assert(mesh.get_dim() == 3);
    assert(mesh.get_vertex_per_face() == 3);
    vertices = mesh.get_vertices().reshape((-1, 3), order="C");
    faces = mesh.get_faces().reshape((-1, 3), order="C");

    vertices, faces, voxels = tetgen(flags, vertices, faces);
    return vertices, faces, voxels;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Tile pattern, tetgen it");
    parser.add_argument("config_file", help="configuration file");
    parser.add_argument("msh_file", help="output msh file");
    args = parser.parse_args();
    return args

def main():
    args = parse_args();

    config = parse_config_file(args.config_file);
    mesh = tile(config);
    vertices, faces, voxels = run_tetgen(config, mesh);
    save_mesh_raw(args.msh_file, vertices, faces, voxels);

if __name__ == "__main__":
    main();


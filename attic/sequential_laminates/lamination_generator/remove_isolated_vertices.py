#!/usr/bin/env python
import numpy as np
import argparse

from mesh_io import load_mesh, form_mesh, save_mesh_raw

def remove_isolated_vertices(vertices, faces, voxels=np.array([], dtype=int)):
    if len(faces) == 0 and len(voxels) == 0:
        return np.array([]), faces, voxels;

    v_set = set(faces.ravel()) | set(voxels.ravel());
    v_set = np.array(list(v_set));
    v_map = np.zeros(len(vertices), dtype=int) - 1;
    v_map[v_set] = np.arange(len(v_set), dtype=int);

    vertices = vertices[v_set];
    faces = v_map[faces];
    voxels = v_map[voxels];
    return vertices, faces, voxels;

def parse_args():
    parser = argparse.ArgumentParser(description="Remove isolated vertices");
    parser.add_argument("in_mesh", help="input mesh");
    parser.add_argument("out_mesh", help="output mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.in_mesh);

    vertices = mesh.vertices;
    faces = mesh.faces;
    voxels = mesh.voxels;

    vertices, faces, voxels = remove_isolated_vertices(vertices, faces, voxels);

    save_mesh_raw(args.out_mesh, vertices, faces, voxels);

if __name__ == "__main__":
    main();

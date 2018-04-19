#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

# Code adapted from quadfoam project, authored by Jeremie Dumas

# System libs
import os
import glob
import json
import struct
import argparse
import subprocess

# Third party libs
import numpy
import pymesh

# Local libs
import paths


def simulate(input_mesh, output_mesh, material_file=None, config_file=None):
    simulator = paths.find("Simulate_cli", paths.MESHFEM_DIR)
    if simulator is None:
        raise FileNotFoundError("Simulate_cli")

    b9material = os.path.join(paths.CONFIG_DIR, "materials/b9creator.json")

    if config_file is None:
        uniaxial_load = os.path.join(paths.CONFIG_DIR, "boundary-conditions/uniaxial_2d_y.json")
    else:
        uniaxial_load = config_file

    if material_file is None:
        material_file = b9material

    args = [simulator, '-m', material_file, '-b', uniaxial_load, '-o', output_mesh, input_mesh]
    subprocess.check_call(args)


def deform_mesh(msh_filename, output_mesh, ratio):
    mesh = pymesh.meshio.load_mesh(msh_filename)
    V = mesh.vertices + ratio * mesh.get_vertex_attribute('u')
    pymesh.meshio.save_mesh_raw(output_mesh, V, mesh.faces)


def export_to_bin(msh_filename, bin_filename):
    mesh = pymesh.meshio.load_mesh(msh_filename)
    with open(bin_filename, 'wb') as file:
        # num_vertices num_triangles num_dofs_attrs num_tri_attrs
        file.write(struct.pack('iiii', mesh.num_vertices, mesh.num_faces, 3, 2))
        assert mesh.vertices.dtype == numpy.float64, "Invalid type"
        assert mesh.faces.dtype == numpy.int32, "Invalid type"
        file.write(bytearray(mesh.vertices.flatten(order='F')))
        file.write(bytearray(mesh.faces.flatten(order='F')))
        for v_attr in ('u', 'load', 'Ku'):
            # TODO: Write attr name
            file.write(bytearray(mesh.get_vertex_attribute(v_attr).flatten(order='F')))
        for f_attr in ('strain', 'stress'):
            # TODO: Write attr name
            file.write(bytearray(mesh.get_face_attribute(f_attr).flatten(order='F')))


def export_profile(msh_filename, mesh_id, db_filename):
    mesh = pymesh.meshio.load_mesh(msh_filename)
    x, y = [], []
    ux, uy = [], []
    thres = mesh.bbox[0][1] + 0.9999 * (mesh.bbox[1][1] - mesh.bbox[0][1])
    print(mesh.bbox)
    for v in range(mesh.vertices.shape[0]):
        if mesh.vertices[v, 1] > thres:
            x.append(mesh.vertices[v, 0])
            y.append(mesh.vertices[v, 1])
            ux.append(mesh.get_vertex_attribute('u')[v, 0])
            uy.append(mesh.get_vertex_attribute('u')[v, 1])
    entry = {
        "id": mesh_id,
        "x": x,
        "y": y,
        "ux": ux,
        "uy": uy
    }
    with open(db_filename, 'r') as f:
        try:
            db = json.load(f)
        except json.decoder.JSONDecodeError:
            db = []
    db.append(entry)
    with open(db_filename, 'w') as f:
        f.write(json.dumps(db, sort_keys=True, indent=4, separators=(',', ': ')))


def process_folder(input_folder, msh_folder, deformed_folder, ratio):
    # Ensure output folders exist
    if not os.path.exists(msh_folder):
        os.makedirs(msh_folder)
    if deformed_folder is not None:
        if not os.path.exists(deformed_folder):
            os.makedirs(deformed_folder)

    # Process all meshes in the input folder
    all_meshes = glob.glob(os.path.join(input_folder, "*.obj"))
    for mesh in all_meshes:
        msh_filename = os.path.splitext(os.path.basename(mesh))[0] + ".msh"
        msh_filename = os.path.join(msh_folder, msh_filename)
        try:
            pass
            # simulate(mesh, msh_filename)
        except subprocess.CalledProcessError:
            print("Could not simulate mesh: " + mesh)
            continue
        if deformed_folder is not None:
            deformed_filename = os.path.splitext(os.path.basename(mesh))[0] + ".ply"
            deformed_filename = os.path.join(deformed_folder, deformed_filename)
            deform_mesh(msh_filename, deformed_filename, ratio)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="input mesh (or folder)")
    parser.add_argument("output", help="output mesh (or folder)")
    parser.add_argument("-m", "--material", type=str, default=None, help="material file (json)")
    parser.add_argument("-e", "--export", type=str, help="filename to export to (for cpp interop)")
    parser.add_argument("-f", "--folder", default=False, action='store_true',
                        help="explicitly enable folder mode (batch processing)")
    parser.add_argument("-d", "--deformed", type=str, help="output deformed mesh")
    parser.add_argument("-r", "--ratio", type=float, default=0.1, help="weight of the deformation field")
    parser.add_argument("-c", "--config", type=str, default=None, help="config file")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.folder or os.path.isdir(args.input):
        process_folder(args.input, args.output, args.deformed, args.ratio)
    else:
        simulate(args.input, args.output, args.material, args.config)
        if args.deformed is not None:
            deform_mesh(args.output, args.deformed, args.ratio)
        if args.export is not None:
            export_to_bin(args.output, args.export)
        # export_profile(args.output, mesh_id=args.input, db_filename="foo.json")


if __name__ == "__main__":
    main()

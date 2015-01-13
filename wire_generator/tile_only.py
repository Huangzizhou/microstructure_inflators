#!/usr/bin/env python

import argparse
import numpy as np
import os.path

from tile import parse_config_file, load_mesh, save_mesh, load_wire, load_parameters, load_wires
from inflator.PyWiresTiler import PyWiresTiler
from wire_io.WireWriter import WireWriter
from core.WireNetwork import WireNetwork

import PyMesh
import PyWires

def save_wire(wires, filename):
    writer = WireWriter(filename);
    writer.write(wires.vertices, wires.edges);

def tile_only(config):
    if "guide_mesh" in config:
        guide_mesh = load_mesh(config["guide_mesh"]);
        if "wire_list_file" in config:
            networks = load_wires(str(config["wire_list_file"]));
            out_wires = tile_with_mixed_patterns(networks, guide_mesh,
                    config.get("dof_type", "isotropic"),
                    config.get("thickness_type", "vertex"));
        else:
            network = load_wire(str(config["wire_network"]));
            parameters = load_parameters(network, config);
            out_wires = tile_with_guide_mesh(network, guide_mesh, parameters);
    else:
        network = load_wire(str(config["wire_network"]));
        parameters = load_parameters(network, config);
        out_wires = tile_with_guide_box(network,
                config["bbox_min"], config["bbox_max"], config["repeats"],
                parameters);
        guide_mesh = generate_box_guide_mesh(
                out_wires.dim,
                np.array(config["bbox_min"]),
                np.array(config["bbox_max"]),
                np.array(config["repeats"], dtype=int));

    offset = out_wires.get_attribute("vertex_offset");
    out_wires.vertices = out_wires.vertices + offset;

    if config.get("trim", False):
        out_wires.trim();
    return out_wires, guide_mesh;

def tile_with_guide_box(unit_pattern, bbox_min, bbox_max, repetitions,
        parameters):
    bbox_min = np.array(bbox_min);
    bbox_max = np.array(bbox_max);
    repetitions = np.array(repetitions, dtype=int);

    tiler = PyWiresTiler();
    tiler.set_single_cell_from_wire_network(unit_pattern);
    tiler.tile(bbox_min, bbox_max, repetitions, parameters);

    tiled_network = tiler.wire_network;

    tiled_bbox_min, tiled_bbox_max = tiled_network.bbox;
    tiled_network.translate(bbox_min - tiled_bbox_min);
    return tiled_network;

def tile_with_guide_mesh(unit_pattern, mesh, parameters):
    tiler = PyWiresTiler();
    tiler.set_single_cell_from_wire_network(unit_pattern);
    tiler.tile_with_hex_mesh(mesh, parameters);
    tiled_network = tiler.wire_network;
    return tiled_network;

def tile_with_mixed_patterns(patterns, mesh, dof_type, thickness_type):
    patterns = [pattern.raw_wires for pattern in patterns];
    tiler = PyWires.WireTiler(patterns[0]);
    tiled_raw_wires = tiler.tile_with_mixed_patterns(
            patterns, mesh,
            thickness_type == "vertex",
            dof_type == "isotropic");

    tiled_wire_network = WireNetwork();
    tiled_wire_network.load_from_raw(tiled_raw_wires);
    return tiled_wire_network;

def generate_box_guide_mesh(dim, bbox_min, bbox_max, num_cells):
    bbox_min = bbox_min[:dim];
    bbox_max = bbox_max[:dim];
    num_cells = num_cells[:dim];
    if dim == 2:
        return generate_box_guide_mesh_2D(bbox_min, bbox_max, num_cells);
    elif dim == 3:
        return generate_box_guide_mesh_3D(bbox_min, bbox_max, num_cells);
    else:
        raise NotImplementedError("Dim={} is not supported".format(dim));

def generate_box_guide_mesh_2D(bbox_min, bbox_max, num_cells):
    assert(np.all(num_cells > 0));
    x_samples = np.linspace(bbox_min[0], bbox_max[0], num_cells[0]+1);
    y_samples = np.linspace(bbox_min[1], bbox_max[1], num_cells[1]+1);
    x_coordinates, y_coordinates = np.meshgrid(x_samples, y_samples);
    vertices = np.hstack((
        x_coordinates.reshape((-1, 1), order="C"),
        y_coordinates.reshape((-1, 1), order="C")));

    x_step = 1;
    y_step = num_cells[0]+1;
    quads = [[
         i   *x_step+ j   *y_step,
         i   *x_step+(j+1)*y_step,
        (i+1)*x_step+(j+1)*y_step,
        (i+1)*x_step+ j   *y_step ]
        for i in range(num_cells[0]) for j in range(num_cells[1]) ];
    quads = np.array(quads, dtype=int);
    voxels = np.zeros((0, 8));

    factory = PyMesh.MeshFactory();
    factory.load_data(
            vertices.ravel(order="C"),
            quads.ravel(order="C"),
            voxels.ravel(order="C"), 2, 4, 0);
    return factory.create_shared();

def generate_box_guide_mesh_3D(bbox_min, bbox_max, num_cells):
    assert(np.all(num_cells > 0));
    x_samples = np.linspace(bbox_min[0], bbox_max[0], num_cells[0]+1);
    y_samples = np.linspace(bbox_min[1], bbox_max[1], num_cells[1]+1);
    z_samples = np.linspace(bbox_min[2], bbox_max[2], num_cells[2]+1);
    x_coordinates, y_coordinates, z_coordinates = \
            np.meshgrid(x_samples, y_samples, z_samples);
    vertices = np.hstack((
        x_coordinates.reshape((-1, 1), order="C"),
        y_coordinates.reshape((-1, 1), order="C"),
        z_coordinates.reshape((-1, 1), order="C")));

    x_step = num_cells[2] + 1;
    y_step = (num_cells[2] + 1) * (num_cells[0] +1);
    z_step = 1;
    hexes = [[
         i   *x_step+ j   *y_step+(k+1)*z_step,
         i   *x_step+(j+1)*y_step+(k+1)*z_step,
        (i+1)*x_step+(j+1)*y_step+(k+1)*z_step,
        (i+1)*x_step+ j   *y_step+(k+1)*z_step,
         i   *x_step+ j   *y_step+ k   *z_step,
         i   *x_step+(j+1)*y_step+ k   *z_step,
        (i+1)*x_step+(j+1)*y_step+ k   *z_step,
        (i+1)*x_step+ j   *y_step+ k   *z_step,
        ]
        for i in range(num_cells[0])
        for j in range(num_cells[1])
        for k in range(num_cells[2]) ];
    hexes = np.array(hexes, dtype=int);
    faces = np.zeros((0, 4));
    tmp = [ (i,j,k)
        for i in range(num_cells[0])
        for j in range(num_cells[1])
        for k in range(num_cells[2]) ];

    factory = PyMesh.MeshFactory();
    factory.load_data(
            vertices.ravel(order="C"),
            faces.ravel(order="C"),
            hexes.ravel(order="C"), 3, 4, 8);
    return factory.create_shared();

def parse_args():
    parser = argparse.ArgumentParser(
            description="Tile a wire network (no inflation)");
    parser.add_argument("--guide-mesh",
            help="output guide mesh used as well", action="store_true");
    parser.add_argument("config_file", help="pattern configuration file");
    parser.add_argument("output_wire", help="output wire file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    config = parse_config_file(args.config_file);
    out_wires, guide_mesh = tile_only(config);

    basename, ext = os.path.splitext(args.output_wire);
    if args.guide_mesh:
        mesh_file = basename + "_guide.obj";
        save_mesh(guide_mesh, mesh_file);

    save_wire(out_wires, args.output_wire);

if __name__ == "__main__":
    main();


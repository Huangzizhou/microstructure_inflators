#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path

import PyWireInflator2DSetting
import PyWireInflator2D
import PyMeshSetting
import PyMesh

from WireModifierFactory import WireModifierFactory
from find_file import find_file

def form_mesh(vertices, faces, voxels=np.array([])):
    dim = vertices.shape[1];
    factory = PyMesh.MeshFactory();
    if dim == 3:
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels.ravel(order="C"), 3, 3, 4);
    elif dim == 2:
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels.ravel(order="C"), 2, 3, 4);
    return factory.create();

def save_mesh(filename, mesh, *attributes):
    writer = PyMesh.MeshWriter.create_writer(filename);
    for attr in attributes:
        if not mesh.has_attribute(attr):
            raise KeyError("Attribute {} is not found in mesh".format(attr));
        writer.with_attribute(attr);
    writer.write_mesh(mesh);

def load_mesh(filename):
    factory = PyMesh.MeshFactory();
    factory.load_file(filename);
    factory.drop_zero_dim();
    return factory.create();

def parse_config_file(config_file):
    """ syntax:
    {
        "wire_network": single_cell_wire_network,
        "thickness": float,
        "modifier_file": modifier_file,
        "bbox_min": [min_x, min_y],
        "bbox_max": [max_x, max_y],
        "repeats": [x_reps, y_reps],
        "quad_mesh": quad_mesh,
        "trim": bool,
        "periodic": bool
    }
    """
    config_dir = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);

    def convert_to_abs_path(field_name):
        field = config[field_name];
        if isinstance(field, (unicode, str)):
            config[field_name] = find_file(field, config_dir);

    convert_to_abs_path("wire_network");
    if "quad_mesh" in config:
        convert_to_abs_path("quad_mesh");
    if "modifier_file" in config:
        convert_to_abs_path("modifier_file");
    return config;

def set_uniform_thickness(wire_network, thickness):
    thickness = np.ones(wire_network.num_vertices) * thickness;
    wire_network.attributes.add("vertex_thickness", thickness);

def load_modifiers(modifier_file):
    if modifier_file is None:
        return [];

    modifiers = WireModifierFactory.create_from_file(str(modifier_file));
    return modifiers;

def load_orbits(inflator):
    basename, ext = os.path.splitext(inflator.source_file);
    orbit_file = basename + ".orbit";
    with open(orbit_file, 'r') as fin:
        orbit_config = json.load(fin);
        inflator.vertex_orbits = orbit_config["vertex_orbits"];

def tile(config):
    wire_file = "{}".format(config["wire_network"]);
    assert(os.path.exists(wire_file));
    inflator = PyWireInflator2D.WireInflatorFacade(wire_file);
    inflator.source_file = wire_file;
    load_orbits(inflator);

    modifiers = load_modifiers(config.get("modifier_file", None));
    if "quad_mesh" in config:
        rows, cols, param, scale_factor = tile_quad(config, modifiers, inflator);
    else:
        rows, cols, param, scale_factor = tile_box(config, modifiers, inflator);

    if config.get("trim", False):
        raise NotImplementedError("Trimming is not supported");

    # Convert wire thickness (unit in mm) to relative thickness used by
    # Luigi's code: radius of the wire assuming each cell is of size 1.0.
    inflator.set_dimension(rows, cols);
    for i in range(rows):
        for j in range(cols):
            p = param[i][j];
            if p is not None:
                inflator.set_parameter(i,j,p);

    if config.get("periodic", False):
        inflator.set_max_triangle_area(0.0001);
        inflator.generate_periodic_pattern();
    else:
        inflator.set_max_triangle_area(0.001);
        inflator.generate_tiled_pattern();

    vertices = inflator.get_vertices().reshape((-1, 2), order="C");
    vertices = vertices * scale_factor;
    faces = inflator.get_triangles().reshape((-1, 3), order="C");

    mesh = form_mesh(vertices, faces, np.array([]));

    attr_names = [];
    if config.get("periodic", False):
        num_vertices = mesh.get_num_vertices();
        mesh.add_attribute("vertex_normal");
        normals = mesh.get_attribute("vertex_normal").reshape((num_vertices, -1));
        velocity = inflator.get_boundary_velocity();
        num_params = velocity.shape[1];
        for i in range(num_params):
            attr_name = "normal_velocity_{}".format(i);
            attr_value = normals * velocity[:,i][:,np.newaxis];
            mesh.add_attribute(attr_name);
            mesh.set_attribute(attr_name, attr_value.ravel(order="C"));
            attr_names.append(attr_name);

    save_mesh(str(config["output"]), mesh, *attr_names);

def tile_quad(config, modifiers, inflator):
    cols, rows = config["repeats"];
    quad_mesh = load_mesh(str(config["quad_mesh"]));

    num_faces = quad_mesh.get_num_faces();
    vertices = quad_mesh.get_vertices().reshape((-1, 2), order="C");
    faces = quad_mesh.get_faces().reshape((-1, 4), order="C");

    bbox_min = np.amin(vertices, axis=0);
    bbox_max = np.amax(vertices, axis=0);
    scale_factor = np.divide(bbox_max - bbox_min, [cols, rows]);
    cell_size = np.divide(bbox_max - bbox_min, [cols, rows])

    attr_names = quad_mesh.get_attribute_names();
    attributes = {}
    for name in attr_names:
        val = quad_mesh.get_attribute(name).ravel();
        if len(val) == num_faces:
            attributes[name] = val;

    face_centers = np.mean(vertices[faces], axis=1);

    num_parameters = inflator.get_num_parameters();
    default_parameter = np.zeros(num_parameters);
    thickness_param_mask = np.zeros(num_parameters, dtype=bool)
    for i in range(num_parameters):
        thickness_param_mask[i] = \
                inflator.get_parameter_type(i) == PyWireInflator2D.WireInflatorFacade.THICKNESS;
    default_parameter[thickness_param_mask] = config["thickness"];

    params = [[None for i in range(cols)] for j in range(rows)];

    def index(p):
        return np.floor(np.divide(p - bbox_min, cell_size)).astype(int)[[1, 0]];

    for i,quad in enumerate(faces):
        row, col = index(face_centers[i]);
        attr_map = {str(name):val[i] for name,val in attributes.iteritems()};
        p = np.copy(default_parameter);
        for modifier in modifiers:
            modifier.modify(p, inflator, **attr_map);
        assert(row < rows and col < cols)
        p[thickness_param_mask] *= 0.5 / np.mean(scale_factor);
        params[row][col] = p;

    return rows, cols, params, scale_factor;

def tile_box(config, modifiers, inflator):
    cols, rows = config["repeats"];
    bbox_min = np.array(config["bbox_min"]);
    bbox_max = np.array(config["bbox_max"]);
    scale_factor = np.divide(bbox_max - bbox_min, [cols, rows]);

    num_params = inflator.get_num_parameters();
    default_parameter = np.zeros(num_params);
    thickness_param_mask = np.zeros(num_params, dtype=bool)
    for i in range(num_params):
        thickness_param_mask[i] = \
                inflator.get_parameter_type(i) == PyWireInflator2D.WireInflatorFacade.THICKNESS;
    default_parameter[thickness_param_mask] = config["thickness"];

    for modifier in modifiers:
        modifier.modify(default_parameter, inflator);

    default_parameter[thickness_param_mask] *= 0.5 / np.mean(scale_factor);
    params = [[np.copy(default_parameter) for i in range(cols)] for j in range(rows)];
    return rows, cols, params, scale_factor;

def parse_args():
    parser = argparse.ArgumentParser(description="Tile a given pattern");
    parser.add_argument("--output", "-o", required=True);
    parser.add_argument("config_file", help="pattern configuration file.");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    config = parse_config_file(args.config_file);
    config["output"] = args.output;
    tile(config);

if __name__ == "__main__":
    main();


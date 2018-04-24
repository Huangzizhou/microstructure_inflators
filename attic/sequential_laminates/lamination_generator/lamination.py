#!/usr/bin/env python
import argparse
import json
import numpy as np
from math import radians, cos, sin, sqrt, modf
import os.path
from subprocess import check_call
from mesh_io import load_mesh, save_mesh_raw
from remove_isolated_vertices import remove_isolated_vertices
from collapse_short_edges import collapse_short_edges;

class LaminationGenerator:
    def __init__(self, config_file, width, height):
        self.__width = width;
        self.__height = height;
        self.__boundary_plate_thickness = 0.025;
        self.__parse_config_file(config_file);

    def generate_laminates(self):
        num_laminates = self.__num_layers;
        layers = [];

        if self.single_material:
            # Always add boundary plates to get clear boundaries.
            # Note triangles on boundary plates might be removed as holes.
            layers.append(self.__generate_left_right_plates());
            layers.append(self.__generate_top_bottom_plates());
        else:
            layers.append(self.__generate_boundary_box());

        for i in range(num_laminates):
            layer = self.__generate_layer(i);
            layers.append(layer);
        self.__laminates = np.vstack(layers);

    def triangulate(self):
        basename, ext = os.path.splitext(self.output_file);
        poly_file = basename + ".poly";
        num_vertices = len(self.__laminates);
        num_laminates = num_vertices / 4;
        with open(poly_file, 'w') as fout:
            fout.write("{} 2 0 0\n".format(num_vertices));
            for i in range(num_vertices):
                fout.write("{} {} {}\n".format(i,
                    self.__laminates[i, 0],
                    self.__laminates[i, 1]));
            fout.write("{} 0\n".format(4 * num_laminates));
            for i in range(num_laminates):
                fout.write("{} {} {}\n".format(i, 4*i  , 4*i+1));
                fout.write("{} {} {}\n".format(i, 4*i+1, 4*i+2));
                fout.write("{} {} {}\n".format(i, 4*i+2, 4*i+3));
                fout.write("{} {} {}\n".format(i, 4*i+3, 4*i  ));
            fout.write("0\n");

        command = "triangle -qpa{:.12f} {}".format(self.max_triangle_area, poly_file);
        check_call(command.split());

        node_file = basename + ".1.node";
        out_file = basename + ext;
        command = "meshconvert.py {} {}".format(node_file, out_file);
        check_call(command.split());

        if self.single_material:
            self.__poke_holes(out_file);
        else:
            self.__trim_to_boundary(out_file);

    def __poke_holes(self, out_file):
        self.__mesh = load_mesh(out_file);
        assert(self.__mesh.vertex_per_face == 3);
        vertices = self.__mesh.vertices;
        faces = self.__mesh.faces;
        face_centers = (
                vertices[faces[:, 0]] +
                vertices[faces[:, 1]] +
                vertices[faces[:, 2]]) / 3.0;
        is_not_hole = np.zeros(len(faces), dtype=bool);
        for i in range(self.__num_layers):
            is_not_hole_i = [not self.__is_hole(p, i) for p in face_centers];
            is_not_hole = np.logical_or(is_not_hole, is_not_hole_i);

        if self.with_top_bottom_plates:
            on_top_bottom_plates = [self.__on_top_bottom_plates(p) for p in face_centers];
            is_not_hole = np.logical_or(is_not_hole, on_top_bottom_plates);
        if self.with_left_right_plates:
            on_left_right_plates = [self.__on_left_right_plates(p) for p in face_centers];
            is_not_hole = np.logical_or(is_not_hole, on_left_right_plates);

        faces = faces[is_not_hole];
        vertices, faces = collapse_short_edges(vertices, faces, 0.0005);
        vertices, faces, voxel = remove_isolated_vertices(vertices, faces);
        save_mesh_raw(out_file, vertices, faces);

    def __trim_to_boundary(self, out_file):
        self.__mesh = load_mesh(out_file);
        assert(self.__mesh.vertex_per_face == 3);
        vertices = self.__mesh.vertices;
        faces = self.__mesh.faces;
        face_centers = (
                vertices[faces[:, 0]] +
                vertices[faces[:, 1]] +
                vertices[faces[:, 2]]) / 3.0;

        half_width = self.width * 0.5;
        half_height = self.height * 0.5;
        is_inside = np.ones(len(faces), dtype=bool);
        is_inside = np.logical_and(is_inside, face_centers[:,0] <  half_width);
        is_inside = np.logical_and(is_inside, face_centers[:,0] > -half_width);
        is_inside = np.logical_and(is_inside, face_centers[:,1] <  half_height);
        is_inside = np.logical_and(is_inside, face_centers[:,1] > -half_height);

        faces = faces[is_inside];
        vertices, faces = collapse_short_edges(vertices, faces, 0.0005);
        vertices, faces, voxel = remove_isolated_vertices(vertices, faces);
        save_mesh_raw(out_file, vertices, faces);

    def __generate_layer(self, index):
        scale = self.scales[index];
        num_laminates = int(max(self.width / scale, self.height / scale));
        laminates = [];
        for i in range(-num_laminates, num_laminates):
            laminate = self.__generate_laminate(index, i);
            laminates.append(laminate);

        return np.vstack(laminates);

    def __generate_laminate(self, layer_idx, lam_idx):
        center = self.centers[layer_idx];
        scale = self.scales[layer_idx];
        ratio = self.ratios[layer_idx];
        thickness = scale * ratio;
        length = 2*sqrt(self.width**2 + self.height**2);

        angle = self.angles[layer_idx];
        direction = self.__rotate(angle, np.array([1, 0]));
        offset = center + direction * lam_idx * scale;

        x = thickness * 0.5;
        y = length * 0.5;
        corners = np.array([[-x, -y], [x, -y], [x, y], [-x, y]]).T;
        corners = self.__rotate(angle, corners).T;
        corners = corners + offset;
        return np.vstack(corners);

    def __generate_top_bottom_plates(self):
        thickness = self.__boundary_plate_thickness;
        x = self.width / 2.0;
        y = self.height / 2.0;
        top = np.array([
            [-x, y-thickness],
            [ x, y-thickness],
            [ x, y+thickness],
            [-x, y+thickness]]);

        y = -self.height / 2.0;
        bottom = np.array([
            [-x, y-thickness],
            [ x, y-thickness],
            [ x, y+thickness],
            [-x, y+thickness]]);
        return np.vstack((top, bottom));

    def __generate_left_right_plates(self):
        thickness = self.__boundary_plate_thickness;
        x = self.width / 2.0;
        y = self.height / 2.0 - 0.1;
        right = np.array([
            [x, -y],
            [x+thickness, -y],
            [x+thickness,  y],
            [x,  y]]);

        x = -self.width / 2.0;
        left = np.array([
            [x-thickness, -y],
            [x, -y],
            [x,  y],
            [x-thickness,  y]]);
        return np.vstack((left, right));

    def __generate_boundary_box(self):
        x = self.width / 2.0;
        y = self.height / 2.0;
        box = np.array([
            [-x,-y],
            [ x,-y],
            [ x, y],
            [-x, y]]);
        return box;

    def __parse_config_file(self, config_file):
        """ syntax:
        {
            "max_triangle_area": #,
            "output": output_filename,
            "with_top_bottom_plates": bool (optional),
            "with_left_right_plates": bool (optional),
            "single_material": bool (optional),
            "laminates": [
                {
                    "center": [x,y],
                    "rotation": angle in degrees,
                    "dimensions": [scale, material ratio]
                }
            ]
        }
        """
        basename, ext = os.path.splitext(config_file);
        path, name = os.path.split(basename);
        with open(config_file, 'r') as fin:
            lamination_spec = json.load(fin);

        self.__num_layers = len(lamination_spec["laminates"]);
        self.__max_triangle_area = lamination_spec["max_triangle_area"];
        self.__output_file = lamination_spec["output"];
        self.__centers = [];
        self.__angles = [];
        self.__scales = [];
        self.__ratios = [];

        if not os.path.isabs(self.__output_file):
            self.__output_file = os.path.join(path, self.__output_file);

        for spec in lamination_spec["laminates"]:
            self.__centers.append(spec["center"]);
            self.__angles.append(radians(spec["rotation"]));
            self.__scales.append(spec["dimensions"][0]);
            self.__ratios.append(spec["dimensions"][1]);

        self.__with_top_bottom_plates = lamination_spec.get("with_top_bottom_plates", False);
        self.__with_left_right_plates = lamination_spec.get("with_left_right_plates", False);
        self.__single_material = lamination_spec.get("single_material", True);

    def __is_hole(self, point, layer_idx):
        if point[0] < -self.width/2.0 or point[0] > self.width/2.0:
            return True;
        if point[1] < -self.height/2.0 or point[1] > self.height/2.0:
            return True;

        center = self.centers[layer_idx];
        scale = self.scales[layer_idx];
        ratio = self.ratios[layer_idx];
        half_ratio = ratio / 2.0;

        angle = self.angles[layer_idx];
        direction = self.__rotate(angle, np.array([1, 0]));

        proj = (point - center).dot(direction);
        frac = modf(proj / scale)[0];
        if frac < 0.0:
            frac = 1 + frac;
        if frac > half_ratio and frac < 1.0 - half_ratio:
            return True;
        else:
            return False;

    def __on_top_bottom_plates(self, point):
        plate_thickness = self.__boundary_plate_thickness;
        half_width = self.width / 2.0;
        half_height = self.height / 2.0;
        if point[0] < -half_width or point[0] > half_width:
            return False;

        on_top = (point[1] <= half_height + plate_thickness) and\
                (point[1] >= half_height - plate_thickness);
        on_bottom = (point[1] <= -half_height + plate_thickness) and\
                (point[1] >= -half_height - plate_thickness);
        return on_top or on_bottom;

    def __on_left_right_plates(self, point):
        plate_thickness = self.__boundary_plate_thickness;
        half_width = self.width / 2.0;
        half_height = self.height / 2.0;
        if point[1] < -half_height or point[1] > half_height:
            return False;

        on_left = (point[0] <= -half_width + plate_thickness) and\
                (point[0] >= -half_width - plate_thickness);
        on_right = (point[0] <= half_width + plate_thickness) and\
                (point[0] >= half_width - plate_thickness);

        return on_left or on_right;


    def __rotate(self, angle, vector):
        return np.array([
            cos(angle) * vector[0] - sin(angle) * vector[1],
            sin(angle) * vector[0] + cos(angle) * vector[1] ]);

    @property
    def centers(self):
        return self.__centers;

    @property
    def angles(self):
        return self.__angles;

    @property
    def scales(self):
        return self.__scales;

    @property
    def ratios(self):
        return self.__ratios;

    @property
    def width(self):
        return self.__width;

    @property
    def height(self):
        return self.__height;

    @property
    def output_file(self):
        return self.__output_file;

    @property
    def max_triangle_area(self):
        return self.__max_triangle_area;

    @property
    def with_left_right_plates(self):
        return self.__with_left_right_plates;

    @property
    def with_top_bottom_plates(self):
        return self.__with_top_bottom_plates;

    @property
    def single_material(self):
        return self.__single_material;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Generate rank-p lamination triangular mesh");
    parser.add_argument("config_file",\
            help="Configure specifying lamination parameters");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    lam_gen = LaminationGenerator(args.config_file, 5, 5);
    lam_gen.generate_laminates();
    lam_gen.triangulate();

if __name__ == "__main__":
    main();

#!/usr/bin/env python

import argparse
import json
import numpy as np
import os.path
from subprocess import check_call
from Holes import Holes 

def parse_config(config_file):
    dir_name = os.path.dirname(config_file);
    with open(config_file, 'r') as fin:
        config = json.load(fin);

    if not os.path.isabs(config["output"]):
        config["output"] = os.path.join(dir_name, config["output"]);
    return config;

def save_poly(output_file, vertices, edges, hole_centers):
    basename, ext = os.path.splitext(output_file);
    poly_file = basename + ".poly";
    with open(poly_file, 'w') as fout:
        fout.write("{} 2 0 0\n".format(len(vertices)));
        for i,v in enumerate(vertices):
            fout.write("{i} {v[0]} {v[1]}\n".format(i=i, v=v));
        fout.write("{} 0\n".format(len(edges)));
        for i,e in enumerate(edges):
            fout.write("{i} {e[0]:d} {e[1]:d}\n".format(i=i, e=e));
        fout.write("{}\n".format(len(hole_centers)));
        for i,h in enumerate(hole_centers):
            fout.write("{i} {h[0]} {h[1]}\n".format(i=i, h=h));
    return poly_file;

def triangulate(poly_file, output_file):
    basename, ext = os.path.splitext(poly_file);
    command = "triangle -pqa0.001 {}".format(poly_file);
    check_call(command.split());

    node_file = basename + ".1.node";
    command = "meshconvert.py {} {}".format(node_file, output_file);
    check_call(command.split());

def generate_holey_rectangle(config):
    """ syntax:
    {
        "lp"       : norm type,
        "radius"   : hole radius,
        "width"    : rectangle width,
        "height"   : rectangle height,
        "grid_size": [num_holes_x, num_holes_y],
        "output"   : output.obj
    }
    """
    p = config["lp"];
    radius = config["radius"];
    width = config["width"];
    height = config["height"];
    grid_size = config["grid_size"];

    dx = width / grid_size[0];
    dy = height / grid_size[1];
    hole_centers = np.array([[i*dx, j*dy]
            for i in range(grid_size[0])
            for j in range(grid_size[1])]);
    hole_centers = hole_centers + [0.5*dx, 0.5*dy];

    holes = Holes(width, height, radius, p);
    holes.drill_holes(hole_centers);

    vertices = np.array([
        [0.0, 0.0],
        [width, 0.0],
        [width, height],
        [0.0, height] ]);
    edges = np.array([
        [0,1],
        [1,2],
        [2,3],
        [3,0] ]);

    vertices = np.vstack((vertices, holes.vertices));
    edges = np.vstack((edges, holes.edges + 4));
    poly_file = save_poly(config["output"], vertices, edges, hole_centers);
    triangulate(poly_file, config["output"]);

def parse_args():
    parser = argparse.ArgumentParser(
            description="Generate rectangle mesh with lp holes");
    parser.add_argument("config_file", help="configuration file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    config = parse_config(args.config_file);
    generate_holey_rectangle(config);

if __name__ == "__main__":
    main();


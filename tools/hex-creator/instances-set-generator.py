#!/usr/bin/env python
import math
import os
import sys
import numpy as np
from subprocess import call


def max_thickness(n, s):
    min_void = 1e-2
    return (s - min_void * (n - 1)) / (2 * n - 1)


def hex_pillars_generator():
    num_pillars_values = range(10, 110, 10)
    triangle_side_values = np.arange(1.0, 1.8, 0.1)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    for num_pillars in num_pillars_values:
        for triangle_side in triangle_side_values:
            max = max_thickness(num_pillars, triangle_side)
            min = max / 2

            thickness_values = np.linspace(min, max, 20)

            for thickness in thickness_values:
                cwd = os.getcwd()
                name = folder_path + '/hexagon-pillars-n{}-s{}-t{}'.format(num_pillars, triangle_side, thickness)
                wire_name = name + '.wire'
                mesh_name = name + '.msh'

                if os.path.isfile(mesh_name):
                    continue

                cmd = [cwd + '/hex-creator.py', str(num_pillars), str(triangle_side), str(thickness), wire_name,
                       mesh_name]
                call(cmd)


def six_leaf_clover_generator():
    center_thickness_parameter_values = np.linspace(0.75, 0.99, 10)
    width_parameter_values = np.linspace(0.5, 0.98, 10)
    tiling_thickness_parameter_values = np.linspace(0.75, 0.99, 10)

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    for center_thickness_parameter in center_thickness_parameter_values:
        for width_parameter in width_parameter_values:
            for tiling_thickness_parameter in tiling_thickness_parameter_values:
                cwd = os.getcwd()
                name = folder_path + '/six-leaf-clover-c{}-w{}-t{}'.format(round(center_thickness_parameter, 3),
                                                                           round(width_parameter, 3),
                                                                           round(tiling_thickness_parameter, 3))
                wire_name = name + '.wire'
                mesh_name = name + '.msh'

                if os.path.isfile(mesh_name):
                    continue

                cmd = [cwd + '/six-leaf-clover-creatorls.py', str(center_thickness_parameter), str(width_parameter),
                       str(tiling_thickness_parameter), wire_name,
                       mesh_name]
                call(cmd)


def hex_diamond_generator():
    triangle_side_values = np.arange(0.5, 1.8, 0.1)
    diamond_width_values = np.arange(0.1, 1.0, 0.1)

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    for triangle_side in triangle_side_values:
        for thickness in diamond_width_values:
            cwd = os.getcwd()
            name = folder_path + '/hexagon-diamond-s{}-t{}'.format(triangle_side, thickness)
            wire_name = name + '.wire'
            mesh_name = name + '.msh'

            if os.path.isfile(mesh_name):
                continue

            cmd = [cwd + '/hexa-diamond-creator.py', str(triangle_side), str(thickness), wire_name,
                   mesh_name]
            call(cmd)


if len(sys.argv) != 3:
    print "usage: ./instances-set-generator.py <generator> <output folder>"
    print "example: ./instances-set-generator.py six-leaf-clover instances"
    sys.exit(-1)

generator = sys.argv[1]
out_path = sys.argv[2]
folder_path = os.getcwd() + '/' + out_path

generator_dictionary = {"hexagon-pillars": hex_pillars_generator,
                        "six-leaf-clover": six_leaf_clover_generator,
                        "hexagon-diamond": hex_diamond_generator}

generator_dictionary[generator]()

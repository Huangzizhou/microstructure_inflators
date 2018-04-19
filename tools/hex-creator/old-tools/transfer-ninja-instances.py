#!/usr/bin/env python
import os
import subprocess
import numpy as np
import sys
import shutil

import hexlib


def print_experiment_info(p1, p2, p3, p4, p5):
    print "\nExperimenting with parameters: \n" \
          "p1: " + str(p1) + "\n" \
                             "p2: " + str(p2) + "\n" \
                                                "p3: " + str(p3) + "\n" \
                                                                   "p4: " + str(p4) + "\n" \
                                                                                "p5: " + str(p5)


def copy_ninja_experiment(triangle_side_factor, num_pillars, chirality_factor, thickness_ratio, ninja_factor):
    print_experiment_info(triangle_side_factor, num_pillars, chirality_factor, thickness_ratio, ninja_factor)
    name = folder_path + '/ninja_p1-{}_p2-{}_p3-{}_p4-{}_p5-{}'.format(round(triangle_side_factor, 3),
                                                                            round(num_pillars, 3),
                                                                            round(chirality_factor, 3),
                                                                            round(thickness_ratio, 3),
                                                                            round(ninja_factor, 3))

    new_name = new_folder_path + '/ninja_p1-{}_p2-{}_p3-{}_p4-{}_p5-{}'.format(round(triangle_side_factor, 3),
                                                                            round(num_pillars, 3),
                                                                            round(chirality_factor, 3),
                                                                            round(thickness_ratio, 3),
                                                                            round(ninja_factor, 3))

    mesh_name = name + '.msh'
    new_mesh_name = new_name + '.msh'

    wire_name = name + '.wire'
    new_wire_name = new_name + '.wire'

    param_name = name + '.wire'
    new_param_name = new_name + '.wire'

    try:
        # free experiment
        shutil.copyfile(mesh_name, new_mesh_name)
    except:
        print "File " + mesh_name + " does NOT exist"
        pass



if len(sys.argv) != 3:
    print "usage: ./ninja-sweep.py <input folder> <output folder>"
    print "example: ./ninja-sweep.py ninja-instances new-ninja-instances"
    sys.exit(-1)

folder_path = sys.argv[1]
new_folder_path = sys.argv[2]

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

chirality_factor_values = [0.78, 0.8, 0.83, 0.85, 0.88, 0.9, 0.93, 0.95]
for chirality_factor in chirality_factor_values:
    ninja_factor_values = [0.7, 0.8, 0.9, 0.95, 1.0]
    for ninja_factor in ninja_factor_values:
        num_pillar_values = [5, 10]
        for index, number_pillars in enumerate(num_pillar_values):
            triangle_side_values = np.arange(0.9, 0.99, 0.02)
            for index, triangle_side_factor in enumerate(triangle_side_values):
                thickness_values = np.arange(0.6, 0.99, 0.04)
                for thickness_factor in thickness_values:
                    copy_ninja_experiment(triangle_side_factor, number_pillars, chirality_factor, thickness_factor, ninja_factor)
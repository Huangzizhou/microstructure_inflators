#!/usr/bin/env python
import math
import os
import sys
import numpy as np
from subprocess import call


def max_thickness(n, s):
    min_void = 1e-2
    return 0.95 * (s - min_void * (n - 1)) / n

if len(sys.argv) != 2:
    print "usage: ./comparison-pillar-types.py <output folder>"
    print "example: ./comparison-pillar-types.py instances"
    sys.exit(-1)

out_path = sys.argv[1]
folder_path = os.getcwd() + '/' + out_path

num_pillars_values = range(5, 25, 5)
triangle_side_values = np.arange(0.8, 1.7, 0.1)
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
for num_pillars in num_pillars_values:
    for triangle_side in triangle_side_values:
        max = max_thickness(num_pillars, triangle_side)
        min = max / 2

        thickness_values = np.linspace(min, max, 20)
        print thickness_values

        for thickness in thickness_values:
            cwd = os.getcwd()
            name_pillars = folder_path + '/hexagon-pillars-n{}-s{}-t{}'.format(num_pillars, triangle_side, thickness)
            name_diamonds = folder_path + '/hexagon-diamonds-n{}-s{}-t{}'.format(num_pillars, triangle_side, thickness)
            name_spears = folder_path + '/hexagon-spears-n{}-s{}-t{}'.format(num_pillars, triangle_side, thickness)

            wire_name = name_pillars + '.wire'
            mesh_name = name_pillars + '.msh'
            if os.path.isfile(mesh_name):
                continue

            cmd = [cwd + '/hexa-many-pillars-creator.py', str(triangle_side), str(num_pillars), str(thickness),
                   wire_name, mesh_name]
            call(cmd)

            wire_name = name_diamonds + '.wire'
            mesh_name = name_diamonds + '.msh'
            if os.path.isfile(mesh_name):
                continue

            cmd = [cwd + '/hexa-many-diamonds-creator.py', str(triangle_side), str(num_pillars), str(thickness),
                   wire_name, mesh_name]
            call(cmd)

            wire_name = name_spears + '.wire'
            mesh_name = name_spears + '.msh'
            if os.path.isfile(mesh_name):
                continue
            cmd = [cwd + '/hexa-many-spears-creator.py', str(triangle_side), str(num_pillars), str(thickness),
                   wire_name, mesh_name]
            call(cmd)

if folder_path.endswith('/'):
    folder_path = folder_path[:-1]
table = folder_path + "-table.txt"
cmd = [cwd + '/run-homogenization.py', folder_path, table]
call(cmd)


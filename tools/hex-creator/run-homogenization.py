#!/usr/bin/env python
import math
import os
import re
import sys
import numpy as np
from subprocess import call

def find_between(s, first, last):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

if len(sys.argv) < 3:
    print "usage: ./run-homogenization.py <input folder> <output table> <material>"
    print "example: ./run-homogenization.py instances table.txt"
    sys.exit(-1)

input_folder = sys.argv[1]
output_table = sys.argv[2]
if len(sys.argv) == 4:
    material = sys.argv[3]
else:
    material = '../../materials/B9Creator.material'


cwd = os.getcwd()
deformed_cells_executable_path = cwd + '/../../../MeshFEM/DeformedCells_cli'

for filename in os.listdir(input_folder):
    if not filename.endswith('msh'):
        continue

    print "Computing homogenized elasticity tensor for: ", filename
    if os.path.isfile(output_table) and filename in open(output_table).read():
        print "Already computed!"
        continue

    cmd = [deformed_cells_executable_path, input_folder + '/' + filename, '-m', material,
           '--homogenize', '--jacobian', '1 0.5 0 0.8660', '--transformVersion']

    with open('output_log.txt', 'w') as out_log:
        call(cmd, stdout=out_log)

    table_file = open(output_table, 'a')
    with open('output_log.txt', 'r') as out_log:
        wire_content = out_log.readlines()
        for line in wire_content:
            if line.startswith('Homogenized Moduli:'):
                pattern = re.compile(r'\-?\d+\.?\d*')  # Compile a pattern to capture float values
                floats = [float(i) for i in pattern.findall(line)]  # Convert strings to float

                youngs_module_x = floats[0]
                youngs_module_y = floats[0]
                youngs_module_avg = (youngs_module_x + youngs_module_y) / 2.0

                poisson_ratio = floats[2]
                shear_modulus = floats[3]

                anisotropy = shear_modulus / (youngs_module_avg / (2 * (1 + poisson_ratio)))

                # infer pattern given name
                pattern = "default"
                if "pillars" in filename:
                    pattern = "pillars"
                if "diamond" in filename:
                    pattern = "diamonds"
                if "spears" in filename:
                    pattern = "spears"

                n = find_between(filename, "-n", "-s")
                pattern += "-n" + n

                table_file.write('{} {} {} {} 0 #{}\n'.format(pattern, youngs_module_avg, poisson_ratio, anisotropy, filename))

                break

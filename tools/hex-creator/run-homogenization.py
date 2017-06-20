#!/usr/bin/env python
import math
import os
import re
import sys
import numpy as np
from subprocess import call

if len(sys.argv) != 3:
    print "usage: ./run-homogenization.py <input folder> <output table>"
    print "example: ./run-homogenization.py instances table.txt"
    sys.exit(-1)

input_folder = sys.argv[1]
output_table = sys.argv[2]

cwd = os.getcwd()
deformed_cells_executable_path = cwd + '/../../../MeshFEM/DeformedCells_cli'

for filename in os.listdir(input_folder):
    if not filename.endswith('msh'):
        continue

    print "Computing homogenized elasticity tensor for: ", filename
    if os.path.isfile(output_table) and filename in open(output_table).read():
        print "Already computed!"
        continue

    cmd = [deformed_cells_executable_path, input_folder + '/' + filename, '-m', '../../materials/B9Creator.material',
           '--homogenize', '--jacobian', '1 0.5 0 0.8660', '--transformVersion']

    with open('output_log.txt', 'w') as out_log:
        call(cmd, stdout=out_log)

    table_file = open(output_table, 'a')
    with open('output_log.txt', 'r') as out_log:
        wire_content = out_log.readlines()
        for line in wire_content:
            if line.startswith('Homogenized Moduli:'):
                pattern = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values
                floats = [float(i) for i in pattern.findall(line)]  # Convert strings to float

                youngs_module = floats[0]
                poisson_ratio = floats[2]

                table_file.write('-6 {} {} 1.0 0 #{}\n'.format(youngs_module, poisson_ratio, filename))

                break

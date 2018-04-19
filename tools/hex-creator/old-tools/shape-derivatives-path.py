#!/usr/bin/env python
import sys

import re, os
import uuid
import subprocess

import numpy as np


def run_optimization(iteration, params_path, wire_path, meshing_file_path, target_Nu, target_E):

    cwd = os.getcwd()

    # create parameters string
    with open(params_path, 'r') as params_file:
        # read parameters file
        content = params_file.readlines()
        parameters = []
        for line_index in range(0, len(content)):
            for element in content[line_index].split():
                # print element
                parameters.append(float(element))

    parameters_string = ', '.join(str(param) for param in parameters)

    print parameters_string

    unique_identifier = str(uuid.uuid1())
    folder_path = 'iteration-'+str(iteration) + '_' + unique_identifier
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Create job
    radiusBounds = [0.00000, 1.0]
    offsetBounds = [-0.5, 0.5]
    blendingBounds = [0.00000, 0.1]

    cmd = ['../../pattern_optimization/GenIsosurfaceJob', wire_path, '-r',
           str(radiusBounds[0]) + ',' + str(radiusBounds[1]),
           '-o', str(offsetBounds[0]) + ',' + str(offsetBounds[1]), '-b',
           str(blendingBounds[0]) + ',' + str(blendingBounds[1]),
           '-e', str(target_E) + ',' + str(target_Nu), '--symmetry', '2D_doubly_periodic', '--initialParams',
           parameters_string]

    job_path = folder_path + '/custom.job'
    with open(job_path, 'w') as out_job:
        subprocess.call(cmd, stdout=out_job)
        print cmd

    # Run optimization
    cmd = ['../../worst_case_stress/WCSOptimization_cli', '-p', wire_path, '-m', '../../materials/Russia.material',
           job_path, '-M', meshing_file_path, '--symmetry', 'doubly_periodic', '--vertexThickness',
           '--WCSWeight', str(0),
           '--JSWeight', str(100.0),
           '--solver', 'levenberg_marquardt',
           '--nIters', str(1),
           #'--solver', 'slsqp',
           '-o', folder_path + '/it',
           '--deformedCell', '1 0.5 0 0.8660']

    out_path = folder_path + '/out.log'
    with open(out_path, 'w') as out_log:
        subprocess.call(cmd, stdout=out_log)
        print cmd

    return out_path


def parse_info(log_path, shape_name, target_Nu, target_E, out_path):
    shape_derivative_norm = -1.0

    with open(log_path, 'r') as log_file:
        for line in log_file:

            m = re.search(r'SHAPE DERIVATIVE NORM: ([0-9]*.[0-9]*)', line)
            if m:
                shape_derivative_norm = float(m.group(1))
                break
            else:
                shape_derivative_norm = -1.0


    with open(out_path, "a") as out_file:
        out_file.write('{} {} {} #{}\n'.format(target_Nu, target_E, shape_derivative_norm, shape_name))


if len(sys.argv) != 8:
    print "usage: ./shape-derivatives-path <initial poisson> <initial youngs module> <final poisson> <final youngs module> <input table> <input folder> <output table>"
    print "example: /shape-derivatives-path 0.0 0.0 -0.5 0.5 ninja-results.txt ninja-results output-table.txt"
    sys.exit(-1)

initial_nu = float(sys.argv[1])
initial_E = float(sys.argv[2])
final_nu = float(sys.argv[3])
final_E = float(sys.argv[4])
table_path = sys.argv[5]
input_folder = sys.argv[6]
out_path = sys.argv[7]

# parsing information in Lookup tables and adding them to the chart
i = 0
nu_values = []
E_values = []
files = []
anysotropy = []
p = []
p1 = []
p2 = []
p3 = []
p4 = []
p5 = []
volfrac = []
table_file = open(table_path)
for line in table_file:
    fields = line.strip().split()
    nu_values.append(float(fields[2]))
    E_values.append(float(fields[1]))
    anysotropy.append(float(fields[3]))
    #files.append(fields[5])
    p.append(fields[0])

    # parse parameters
    m = re.search('volfrac-(.+?)_', fields[5])
    if m:
        volfrac.append(float(m.group(1)))
    else:
        volfrac.append(0.0)

    m = re.search('p1-(.+?)_', fields[5])
    if m:
        p1.append(float(m.group(1)))
    else:
        p1.append(0.0)

    m = re.search('p2-(.+?)_', fields[5])
    if m:
        p2.append(float(m.group(1)))
    else:
        p2.append(0)

    m = re.search('p3-(.+?)_', fields[5])
    if m:
        p3.append(float(m.group(1)))
    else:
        p3.append(0)

    m = re.search('p4-(.+?)_', fields[5])
    if m:
        p4.append(float(m.group(1)))
    else:
        p4.append(0)

    m = re.search('p5-(.+?).msh', fields[5])
    if m:
        p5.append(float(m.group(1)))
    else:
        p5.append(0)

    m = re.search('#(.+?).msh', fields[5])
    if m:
        files.append(m.group(1))
    else:
        files.append(0)

lut = []

for i in range(0,len(nu_values)):
    print i
    lut.append([nu_values[i], E_values[i]])

np_lut = np.array(lut)

initial_point = np.array([initial_nu, initial_E])
final_point = np.array([final_nu, final_E])
direction = final_point - initial_point
steps = 10

for t in range(0, steps+1):
    target = initial_point + 1.0 * t / steps * direction
    print target

    # find E and Nu that are the closest to our current objective
    distToTarget = np_lut - target
    lutDists = np.linalg.norm(distToTarget, axis=1)
    order = np.argsort(lutDists)
    friends = min(1, len(lut))
    shape_name = files[order[0]]
    shape_path = input_folder + '/' + shape_name
    wire_path = shape_path + '.wire'
    params_path = shape_path + '.param'

    log_path = run_optimization(t, params_path, wire_path, 'shape_derivatives_meshing_file.json', final_nu, final_E)

    parse_info(log_path, shape_name, nu_values[order[0]], E_values[order[0]], out_path)


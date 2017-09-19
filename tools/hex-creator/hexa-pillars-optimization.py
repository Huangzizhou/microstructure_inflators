#!/usr/bin/env python
import json
import os
import sys
import subprocess
import uuid

if len(sys.argv) != 8:
    print "usage: ./hexa-pillars-optimization.py <+/-> <p1> <p2> <p3> <p4> <target Poisson ratio> <target Young's module>"
    print "example: ./hexa-pillars-optimization.py + 0.5 6 0.7 0.8 0.25 0.75"
    sys.exit(-1)

type = sys.argv[1]
p1 = sys.argv[2]
p2 = sys.argv[3]
p3 = sys.argv[4]
p4 = sys.argv[5]
target_Nu = sys.argv[6]
target_E = sys.argv[7]

cwd = os.getcwd()

unique_identifier = str(uuid.uuid1())
folder_path = 'temp-' + unique_identifier
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

name = folder_path + '/struct'
wire_name = name + '.wire'
wire_mesh = name + '.msh'
cmd = [cwd + '/hex-creator-igor-parameters.py', str(p1), str(p2), str(p3), str(p4), name + '.wire', name + '.msh']
subprocess.call(cmd)

original_job_path = folder_path + '/original.job'
custom_job_path = folder_path + '/custom.job'
cmd = ['../../pattern_optimization/GenIsosurfaceJob', wire_name, '-e', str(target_E) + ',' + str(target_Nu),
       '--symmetry', '2D_doubly_periodic']

with open(original_job_path, 'w') as out_job:
    ret = subprocess.call(cmd, stdout=out_job)

with open(original_job_path) as original_file:
    job_opts = json.load(original_file)

    if type == "+":
        job_opts['initial_params'] = [p1, p4]
        job_opts['meta_params'] = ['+', int(p2)]
    else:
        job_opts['initial_params'] = [p1, p3, p4]
        job_opts['meta_params'] = ['-', int(p2)]

    job_opts['metaBounds'] = [0.5, 0.95]
    job_opts['custom1Bounds'] = [0.7, 0.98]
    job_opts['custom3Bounds'] = [0.5, 0.95]
    job_opts['custom4Bounds'] = [0.5, 0.95]


    with open(custom_job_path, 'w') as outfile:
        json.dump(job_opts, outfile)

if type == "+":
    parameters_string = str(p1) + ', ' + str(p4)
    meta_parameters_string = '+' + ', ' + str(int(p2))
else:
    parameters_string = str(p1) + ', ' + str(p3) + ', ' + str(p4)
    meta_parameters_string = '-' + ', ' + str(int(p2))

# Run optimization
cmd = ['../../worst_case_stress/WCSOptimization_cli', custom_job_path,
       '--params', parameters_string,
       '--metaParams', meta_parameters_string,
       '-i', 'HexaPillars',
       '-m', '../../materials/Russia.material',
       '--vertexThickness',
       '--WCSWeight', str(0),
       '--JSWeight', str(10.0),
       #'--proximityRegularizationWeight', str(1.0),
       '--solver', 'levenberg_marquardt',
       '-o', folder_path + '/it',
       '--deformedCell', '1 0.5 0 0.8660']
print cmd
out_path = folder_path + '/optimization.log'
with open(out_path, 'w') as out_log:
    ret = subprocess.call(cmd)

# adding results to table, running homogenization on 'it' instances
for filename in os.listdir(folder_path):
    if filename.startswith("it_"):
        os.rename(folder_path + "/" + filename, folder_path + "/" + filename + ".msh")

table_path = folder_path + '/table.txt'
cmd = [cwd + '/run-homogenization.py', folder_path, table_path, "../../materials/Russia.material"]
subprocess.call(cmd)

# plotting results
cmd = [cwd + '/plot_luTable.py', 'all', "../../materials/Russia.material", table_path]
subprocess.call(cmd)


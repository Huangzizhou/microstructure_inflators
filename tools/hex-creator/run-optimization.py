#!/usr/bin/env python
import os
import sys
import subprocess

if len(sys.argv) != 5:
    print "usage: ./run-optimization.py <wire file> <params file> <target E> <target nu>"
    print "example: ./run-optimization.py instance.wire instance.params 10 0.7"
    sys.exit(-1)

wire_path = sys.argv[1]
params_path = sys.argv[2]
target_E = sys.argv[3]
target_Nu = sys.argv[4]

cwd = os.getcwd()

radiusBounds = [0.005, 0.3]
offsetBounds = [-0.4, 0.4]
blendingBounds = [0.001, 0.5]
dim = 2

# Create job
cmd = ['../../pattern_optimization/GenIsosurfaceJob', wire_path, '-r', str(radiusBounds[0]) + ',' + str(radiusBounds[1]),
    '-o', str(offsetBounds[0]) + ',' + str(offsetBounds[1]), '-b', str(blendingBounds[0]) + ',' + str(blendingBounds[1]),
    '-e', target_E + ',' + target_Nu]

job_path = os.path.splitext(wire_path)[0] + '.job'
with open(job_path, 'w') as out_job:
    ret = subprocess.call(cmd, stdout=out_job)

# Run optimization
cmd = ['../../worst_case_stress/WCSOptimization_cli', '-p', wire_path, '-m', '../../materials/B9Creator.material',
       job_path, '-M', 'refined-meshing_opts.json', '--ortho-cell', '--vertexThickness', '-solver', 'slsqp', '-o', 'it']

out_path = os.path.splitext(wire_path)[0] + '.log'
with open(out_path, 'w') as out_log:
    ret = subprocess.call(cmd, stdout=out_log)
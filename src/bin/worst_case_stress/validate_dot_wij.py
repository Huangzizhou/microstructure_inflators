#! /usr/bin/env python
# dotw_ij is supposed to be the Eulerian derivative of w_ij. We validate this on
# a set of sample points over a sequence of meshes using finite differences.
import numpy as np
import subprocess
import pandas as pd

samplePts = [[0.00, 0.55], [0.00, 0.75], [0.00, 0.90], [0.00, 0.95],
             [0.55, 0.55], [0.75, 0.75], [0.90, 0.90], [0.95, 0.95]]
samplePtNames = ['s%i' % i for i in range(len(samplePts))] 

fields = ['w 0', 'w 1', 'w 2', 'p1 wdot 0', 'p1 wdot 1', 'p1 wdot 2']
componentFields = []
for f in fields: componentFields += [f + ' x', f + ' y']

# dict giving dict of sample point values for each field.
data = {cf: {sp: [] for sp in samplePtNames} for cf in componentFields}

sampleString = ";".join(map(lambda p: "{}, {}".format(p[0], p[1]), samplePts))

for i in range(60):
    cmd = ["msh_processor", "lphole_iterates/param_p_%i.msh" % i]
    for f in fields: cmd += ['-e', f]
    cmd += ['-A', '--sample', sampleString, '--reverse', '-Ap']
    out = subprocess.check_output(cmd).strip().split()
    idx = 0
    for cf in componentFields:
        for sp in samplePtNames:
            data[cf][sp].append(float(out[idx]))
            idx = idx + 1
p = pd.Panel(data)
p.to_hdf("data.hdf", "data", mode='w')
# field (including x, y)
# sample point
# iteration

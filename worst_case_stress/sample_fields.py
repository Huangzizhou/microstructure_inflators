#! /usr/bin/env python
# dotw_ij is supposed to be the Eulerian derivative of w_ij. We validate this on
# a set of sample points over a sequence of meshes using finite differences.
import numpy as np
import subprocess
import pandas as pd
import re

nonfloatChars = re.compile(r'[^\d.+\-e]+')
removeNonfloatParts = lambda s: nonfloatChars.sub('', s)

samplePts = [[0.00, 0.55], [0.00, 0.75], [0.00, 0.90], [0.00, 0.95],
             [0.55, 0.55], [0.75, 0.75], [0.90, 0.90], [0.95, 0.95]]
samplePtNames = ['s%i' % i for i in range(len(samplePts))] 

sampleString = ";".join(map(lambda p: "{}, {}".format(p[0], p[1]), samplePts))

componentNames = ['x', 'y', 'xy']

params = ['radius', 'p']
for pnum in [1]:
    param = params[pnum]
    fields = ['w 0', 'w 1', 'w 2']
    numFieldComponents = 3 * [2]

    fields += ['we 0', 'we 1', 'we 2']
    numFieldComponents += 3 * [3]

    fields += ['j']
    numFieldComponents += [1]

    prefix = "p%i " % pnum;
    derivativeFields  = [prefix + s for s in ['wdot 0', 'wdot 1', 'wdot 2']]
    numDerivativeFieldComponents = 3 * [2]

    derivativeFields += [prefix + s for s in ['wedot 0', 'wedot 1', 'wedot 2']]
    numDerivativeFieldComponents += 3 * [3]

    derivativeFields += [prefix + s for s in ['jdot', 'fd jdot', 'corrected fd jdot', 'gradC corrected fd jdot']]
    numDerivativeFieldComponents += 4 * [1]

    # dict giving dict of sample point values for each field.
    data, derivativeData = {}, {}
    for f, nc in zip(fields, numFieldComponents):
        for cf in [f + ' ' + c if nc > 1 else f for c in componentNames[0:nc]]:
            data[cf] = {sp: [] for sp in samplePtNames}
    for f, nc in zip(derivativeFields, numDerivativeFieldComponents):
        for cf in [f + ' ' + c if nc > 1 else f for c in componentNames[0:nc]]:
            derivativeData[cf] = {sp: [] for sp in samplePtNames}

    paramValues = []
    for line in file('lphole_iterates/' + param + '.txt'):
        if (re.match("[0-9]+\t[0-9]", line)):
            paramValues += line.strip().split()[1:2]
    paramValues = map(float, paramValues)
    for i in range(len(paramValues)):
        cmd = ["msh_processor", "lphole_iterates/param_%s_%i.msh" % (param, i)]
        for f in fields: cmd += ['-e', f]
        for f in derivativeFields: cmd += ['-e', f]
        cmd += ['-A', '--sample', sampleString, '--reverse', '-Ap']
        out = subprocess.check_output(cmd).strip().split()
        idx = 0
        for f, nc in zip(fields, numFieldComponents):
            for sp in samplePtNames:
                for c in componentNames[0:nc]:
                    name = f + ' ' + c if nc > 1 else f
                    data[name][sp].append(float(removeNonfloatParts(out[idx])))
                    idx = idx + 1
        for f, nc in zip(derivativeFields, numDerivativeFieldComponents):
            for sp in samplePtNames:
                for c in componentNames[0:nc]:
                    name = f + ' ' + c if nc > 1 else f
                    derivativeData[name][sp].append(float(removeNonfloatParts(out[idx])))
                    idx = idx + 1
        if (idx != len(out)): raise Exception("Unused data")
    store = pd.HDFStore(param + ".hdf")
    datPanel = pd.Panel(data)
    derDatPanel = pd.Panel(derivativeData)

    store['w'] = datPanel.select(lambda label: 'j' not in label)
    store['wdot'] = derDatPanel.select(lambda label: 'j' not in label)
    store['paramValues'] = pd.Series(paramValues)

    store['j'] = datPanel.select(lambda label: 'j' in label)
    store['jdot'] = derDatPanel.select(lambda label: 'j' in label)

    store.close()

#!/bin/env python
# Extracts a component lookup table from a collection of DeformedCells_cli output files.
# This is a matrix of size nthetas x nlambdas where each entry is the value of
# a single tensor component. The component to extract is specified as an argument.
import sys;
from collections import defaultdict;

if len(sys.argv) < 3:
    raise Exception("Usage: extract_component_lut.py component_num resultFile1 [resultFile2 ...]");

component = int(sys.argv[1]);
outfiles = sys.argv[2:];

allLines = [];
for filename in outfiles:
    allLines.extend(list(open(filename)));

# values[theta][lambda] = component_i
values = defaultdict(dict);
for l in allLines:
    row = l.strip().split('\t');
    values[row[0]][row[1]] = row[2:][component];

# Validate table
thetaKeys = sorted(values.keys());
lambdaKeys = sorted(values[thetaKeys[0]].keys());
for k in thetaKeys:
    if sorted(values[k].keys()) != lambdaKeys:
        raise Exception('Invalid table keys');

# Write out table;
for tkey in thetaKeys:
    for lkey in lambdaKeys:
        print values[tkey][lkey],;
    print "";

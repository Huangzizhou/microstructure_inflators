#!/usr/bin/env python
# Parse shape velocity/gradient validation output and print table:
#   paramValue objective/volume shapeVelocityDObj centeredDifffDObj forwardDiffDObj
# assumes the parameter values are evenly spaced.
import sys, re, numpy as np;
dataPath = sys.argv[1]
paramValues = []
objective = []
shapeVelocityDObj = []
for line in file(dataPath):
    if (re.match("[0-9]+\t[0-9]", line)):
        fields = line.strip().split()
        paramValues += [fields[1]]
        objective += [fields[2]]
        shapeVelocityDObj += [fields[3]]
paramValues = np.array(map(float, paramValues))
objective = np.array(map(float, objective))
shapeVelocityDObj = np.array(map(float, shapeVelocityDObj))
centeredDifffDObj = np.gradient(objective, paramValues[1] - paramValues[0])
forwardDiffDObj = np.diff(objective) / np.diff(paramValues)
for i in range(len(paramValues)):
    print "{}\t{}\t{}\t{}\t{}".format(paramValues[i], objective[i], shapeVelocityDObj[i],
                                  centeredDifffDObj[i], "N/A" if i == 0 else forwardDiffDObj[i - 1])

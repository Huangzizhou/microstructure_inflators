#!/usr/bin/env python
# Parse shape velocity/gradient validation output and print table:
#   paramValue objective/volume shapeVelocityDObj centeredDiffDObj forwardDiffDObj wcObjective wcDirectDObj wcAdjointDObj wcCenteredDiffObj wcForwardDiffObj wcDirectKLTerm wcDirectAdvectTerm 
# assumes the parameter values are evenly spaced.
import sys, re, numpy as np;
dataPath = sys.argv[1]

paramValues, objective, shapeVelocityDObj = [], [], []
wcObjective, wcDirectKLTerm, wcDirectAdvectTerm, wcDirectDObj, wcAdjointDObj = [], [], [], [], []
for line in file(dataPath):
    if (re.match("[0-9]+\t[0-9]", line)):
        fields = line.strip().split()
        paramValues += [fields[1]]
        objective += [fields[2]]
        shapeVelocityDObj += [fields[3]]
        wcObjective += [fields[4]]
        wcDirectKLTerm += [fields[5]]
        wcDirectAdvectTerm += [fields[6]]
        wcDirectDObj += [fields[7]]
        wcAdjointDObj += [fields[8]]
paramValues = np.array(map(float, paramValues))
objective = np.array(map(float, objective))
wcObjective = np.array(map(float, wcObjective))

centeredDiffDObj = np.gradient(objective, paramValues[1] - paramValues[0])
forwardDiffDObj = np.diff(objective) / np.diff(paramValues)

wcCenteredDiffDObj = np.gradient(wcObjective, paramValues[1] - paramValues[0])
wcForwardDiffDObj = np.diff(wcObjective) / np.diff(paramValues)

for i in range(len(paramValues)):
    print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(paramValues[i], objective[i], shapeVelocityDObj[i],
                                  centeredDiffDObj[i], "N/A" if i == 0 else forwardDiffDObj[i - 1],
                                  wcObjective[i], wcDirectDObj[i], wcAdjointDObj[i], wcCenteredDiffDObj[i],
                                  "N/A" if i == 0 else wcForwardDiffDObj[i - 1],
                                  wcDirectKLTerm[i], wcDirectAdvectTerm[i])

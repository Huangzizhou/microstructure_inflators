import LookupTable, sys
import numpy as np
import json
from autocover_tools import materialSpaceGrid
from numpy.linalg import norm
import math

lutPath, configPath, outPath = sys.argv[1:]

lut = LookupTable.LUT(lutPath)
config = json.load(file(configPath))
grid = materialSpaceGrid(config)

# compute current covering grid:
currERange = (lut.minE(), lut.maxE());
currNuRange = (lut.minNu(), lut.maxNu());
minIndices = map(int, grid.getPointIndex(currERange[0], currNuRange[0], math.floor))
maxIndices = map(int, grid.getPointIndex(currERange[1], currNuRange[1], math.ceil))


# Compute indices of all patterns
tableUnroundedIndices = np.zeros((lut.size(), 2))
for i, yna in enumerate(lut.youngPoissonAnisotropy):
    tableUnroundedIndices[i] = grid.getPointIndex(yna[0], yna[1], roundingOp=lambda x: x)
    closestPtIndices = map(round, tableUnroundedIndices[i])

# Choose the single closest pattern to each gridpoint
gridPoints = set()

deletePattern = np.ones(lut.size(), dtype=bool)

for EIdx in range(minIndices[0], maxIndices[0] + 1):
    for nuIdx in range(minIndices[1], maxIndices[1] + 1):
        targetE, targetNu = grid.pointAtIndices(EIdx, nuIdx)
        pt = np.array([EIdx, nuIdx])
        lutPointDists = norm(tableUnroundedIndices - pt, axis=1)
        order = np.argsort(lutPointDists)

        # Skip gridpoints too far away from existing points
        if (lutPointDists[order[0]] > 2):
            continue

        deletePattern[order[0]] = False

lut.removeIndices(np.where(deletePattern))

lut.write(outPath)

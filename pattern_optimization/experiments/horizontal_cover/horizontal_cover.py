#!/usr/bin/env python
import json, sys, paths
import numpy as np
from optimization import *
from LookupTable import LUT, extract as extractLUT
import shutil
import glob
import pdb

def deleteFolders(directory, prefix):
    regex = directory + '/' + prefix + '*'
    foldersList = glob.glob(regex)

    print "Deleting following folders:"
    for folder in foldersList:
        print(folder)
        shutil.rmtree(folder)


if (len(sys.argv) != 2):
    print "usage: ./horizontal_cover.py config"
    print "example: ./horizontal_cover.py 2D_autocover_config.json"
    sys.exit(-1)

configPath = sys.argv[1]
config = json.load(file(configPath))

pattern = config['pattern']
radiusBounds = config['radiusBounds']
offsetBounds = config['offsetBounds']
blendingBounds = config['blendingBounds']
translationBounds = config['translationBounds']
climbingSteps = config['climbingSteps']
branches = config['branches']
dim = config['dim']
printable = config['printable']

# Need to implement this to get template job:
# $MICRO_DIR/pattern_optimization/GenIsosurfaceJob $i -e'1,0' -r'0.02,0.2' -o'-0.3,0.3' > init_jobs/$pat.opt;
cwd = os.getcwd() 
cmd = [paths.jobGenerator(config['dim']), paths.pattern(pattern, config['dim']), 
    '-r', str(radiusBounds[0]) + ',' + str(radiusBounds[1]), 
    '-o', str(offsetBounds[0]) + ',' + str(offsetBounds[1]), 
    '-b', str(blendingBounds[0]) + ',' + str(blendingBounds[1]), 
    '-e', str(config['initialE']) + ',' + str(config['initialNu'])]
print(cmd)

outPath = cwd + '/template.job'
with open(outPath, 'w') as outLog:
    ret = subprocess.call(cmd, stdout=outLog)

config['jobTemplate'] = outPath
initialJobPath = outPath
initialJob = json.load(file(initialJobPath))

# Enqueue and execute job
resultsPath = cwd + '/results'
deleteFolders(cwd, 'results')

opt = PatternOptimization(config)
opt.enqueueJob(config['initialE'], config['initialNu'], initialJob['initial_params'])
opt.run(resultsPath)

# Extract table with data achieved
lut = extractLUT(dim, pattern, resultsPath, printableOnly = printable)

it = 1
for curTargetE in np.linspace(config['initialE'], config['targetE'], climbingSteps):
    newTarget = (curTargetE, config['initialNu'], 1)
    print(newTarget)

    # Use look up table to find parameters for a new run, with a bit larger larger targetE
    tableLut = np.zeros((lut.size(), 3))
    for i, yna in enumerate(lut.youngPoissonAnisotropy):
        tableLut[i] = yna

    # find E and Nu that are the closest to our current objective
    distToTarget = tableLut - newTarget
    distToTarget[:,0] *= 1.0/newTarget[0]
    distToTarget[:,1] *= 1.0/newTarget[1]
    lutDists = np.linalg.norm(distToTarget, axis=1)
    order = np.argsort(lutDists)
    
    friends = min(branches, lut.size())
    for i in range(0,friends):
        print "\tUsed to help:"
        print "\t ", (tableLut[order[i]])
        opt.enqueueJob(newTarget[0], newTarget[1], lut.params[order[i]])

    opt.run(resultsPath + '_' + str(it))

    newLut = extractLUT(dim, pattern, resultsPath + '_' + str(it), printableOnly = printable)
    lut.union(newLut)

    it = it + 1
    
lut.filterAnisotropy(0.90, 1.10)
lut.write(cwd + '/climbingLutTable.txt')

# Now, after climbing to our start point. Let's travel to the limit on both sides
NSubdiv = config['targetNSubdiv']
targetNuRange = config['targetNuRange']
targetE = config['targetE']
initialNu = config['initialNu']
initialE = config['initialE']

nuValues = np.linspace(targetNuRange[0], targetNuRange[1], NSubdiv)
orderNu = np.argsort(np.abs(nuValues - initialNu))
interval = nuValues[1] - nuValues[0]
print "interval size: ", interval
resultsPath = cwd + '/results_horizontal'
it = 1
for i in range(0, NSubdiv):
    targetNu = nuValues[orderNu[i]]
    newTarget = (targetE, targetNu, 1) 
    print(newTarget)

    # Use look up table to find parameters for a new run, with a bit larger targetE
    tableLut = np.zeros((lut.size(), 3))
    for i, yna in enumerate(lut.youngPoissonAnisotropy):
        tableLut[i] = yna

    # find E and Nu that are the closest to our current objective
    distToTarget = tableLut - newTarget
    distToTarget[:,0] *= 1.0/newTarget[0]
    distToTarget[:,1] *= 5.0/newTarget[1]

    lutDists = np.linalg.norm(distToTarget, axis=1)
    orderDist = np.argsort(lutDists)
    
    bestDiffNu = tableLut[orderDist[0]][1] - targetNu
    if (abs(bestDiffNu) > 5*interval):
        print "Not worth trying!"
        continue

    for i in range(0,branches):
        print "\tUsed to help:"
        print "\t ", (tableLut[orderDist[i]])
        opt.enqueueJob(newTarget[0], newTarget[1], lut.params[orderDist[i]])

    opt.run(resultsPath + '_' + str(it))

    newLut = extractLUT(dim, pattern, resultsPath + '_' + str(it), printableOnly = True)
    lut.union(newLut)

    it = it + 1

lut.filterAnisotropy(0.90, 1.10)
lut.filterYoung(targetE - 5, targetE + 5)
lut.write(cwd + '/horizontalLutTable.txt')

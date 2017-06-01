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



if (len(sys.argv) != 7):
    print "usage: ./gradient_generator.py config table outPrefix pattern youngsModule poissonRatio"
    print "example: ./gradient_generator.py 2D_autocover_config.json table.txt it 98 20 0.5"
    sys.exit(-1)

configPath = sys.argv[1]
tablePath = sys.argv[2]
outPrefix = sys.argv[3]
pattern = int(sys.argv[4])
targetE = float(sys.argv[5])
targetNu = float(sys.argv[6])

config = json.load(file(configPath))

radiusBounds = config['radiusBounds']
offsetBounds = config['offsetBounds']
blendingBounds = config['blendingBounds']
translationBounds = config['translationBounds']
climbingSteps = config['climbingSteps']
branches = config['branches']
dim = config['dim']
printable = config['printable']
symmetry = config.get('symmetry')
limitedOffset = config.get('limitedOffset')

# set pattern to config file
config['pattern'] = pattern

# create job to run
cwd = os.getcwd() 
cmd = [paths.jobGenerator(config['dim']), paths.pattern(pattern, config['dim']), 
    '-r', str(radiusBounds[0]) + ',' + str(radiusBounds[1]), 
    '-o', str(offsetBounds[0]) + ',' + str(offsetBounds[1]), 
    '-b', str(blendingBounds[0]) + ',' + str(blendingBounds[1]), 
    '-e', str(config['initialE']) + ',' + str(config['initialNu'])]

print(cmd)

if symmetry:
    cmd = cmd + ['--' + str(symmetry)]

if limitedOffset:
    cmd = cmd + ['--limitedOffset']

outPath = cwd + '/template.job'
with open(outPath, 'w') as outLog:
    ret = subprocess.call(cmd, stdout=outLog)

config['jobTemplate'] = outPath
initialJobPath = outPath
initialJob = json.load(file(initialJobPath))

# Extract table
lut = LUT(tablePath)

# Use look up table to find parameters for a new run, with a bit larger larger targetE
tableLut = np.zeros((lut.size(), 3))
for i, yna in enumerate(lut.youngPoissonAnisotropy):
    tableLut[i] = yna
    
newTarget = (targetE, targetNu, 1)

# find E and Nu that are the closest to our current objective
distToTarget = tableLut - newTarget
distToTarget[:,0] *= 1.0/newTarget[0]
distToTarget[:,1] *= 1.0/newTarget[1]
lutDists = np.linalg.norm(distToTarget, axis=1)
order = np.argsort(lutDists)

# Enqueue and execute job
resultsPath = cwd + '/results-gradient' + str(newTarget)
deleteFolders(cwd, 'results-gradient' + str(newTarget))

params = lut.params[order[0]]
closestE = tableLut[order[0]][0]
closestNu = tableLut[order[0]][1]

opt = PatternOptimization(config)
opt.enqueueJob(targetE, targetNu, lut.params[order[0]])
#opt.enqueueJob(closestE, closestNu, lut.params[order[0]])
opt.run(resultsPath)



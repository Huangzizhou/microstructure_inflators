#!/usr/bin/env python
import json, sys, paths
import numpy as np
from optimization import *
from LookupTable import LUT, extract as extractLUT
import shutil
import glob
import pdb

if (len(sys.argv) != 5) and (len(sys.argv) != 6):
    print "usage: ./elasticity_to_mesh.py luTable output youngsModule poissonRatio [inflationGraph]"
    print "example: ./elasticity_to_mesh.py table.txt drawing.msh 20 0.3 inflation.msh"
    sys.exit(-1)

tablePath = sys.argv[1]
outPath = sys.argv[2]
targetE = float(sys.argv[3])
targetNu = float(sys.argv[4])
inflationGraph = ""
dim = 2

if (len(sys.argv) == 6):
    inflationGraph = sys.argv[5]

inflatorPath = paths.inflator()

lut = LUT(tablePath)
tableLut = np.zeros((lut.size(), 3))
for i, yna in enumerate(lut.youngPoissonAnisotropy):
    tableLut[i] = yna

# find E and Nu that are the closest to our current objective
target = (targetE, targetNu, 1)
distToTarget = tableLut - target
distToTarget[:,0] *= 1.0/target[0]
distToTarget[:,1] *= 1.0/target[1]
lutDists = np.linalg.norm(distToTarget, axis=1)
order = np.argsort(lutDists)

for i in range(0,3):
    print "\tClosest:"
    print "\t ", (tableLut[order[i]])

params = lut.params[order[0]]
paramsString = str(params)
paramsStringLen = len(paramsString)
paramsString = paramsString[1:(paramsStringLen-1)].replace('\n', ' ')
print paramsString

# create json file to contain mesh options
meshingOpts = {}
meshingOpts['domainErrorBound'] = 1e-5
meshingOpts['facetAngle'] = 30.0
meshingOpts['facetSize'] = 0.025
meshingOpts['facetDistance'] =  2e-3
meshingOpts['cellSize'] = 0.15
meshingOpts['edgeSize'] = 0.025
meshingOpts['cellRadiusEdgeRatio'] = 2.0
meshingOpts['marchingSquaresGridSize'] = 512
meshingOpts['marchingCubesGridSize'] = 128
meshingOpts['maxArea'] = 5e-4
meshingOpts['featureAngleThreshold'] = 0.7853981633974483
meshingOpts['forceMSGridSize'] = True
meshingOpts['marchingSquaresCoarsening'] = 3

print meshingOpts

meshingOptsPath = "tempMeshingOpts.txt"
with open(meshingOptsPath, 'w') as outfile:
    json.dump(meshingOpts, outfile)

cwd = os.getcwd() 
cmd = [paths.inflator(), '2D_orthotropic', paths.pattern(int(lut.pattern[order[0]]), dim), 
        outPath, '-m',  meshingOptsPath,'--params', paramsString]
print cmd
try:
    subprocess.check_output(cmd)
except Exception, e:
    print "Did not work out perfectly!! Maybe it is not orthotropic!?"
    print str(e)
   
    print "Let's try again with 2D_square!"
    
    if (len(sys.argv) == 6):
        cmd = [paths.inflator(), '2D_square', paths.pattern(int(lut.pattern[order[0]]), dim), 
                outPath, '-m',  meshingOptsPath, '-D', inflationGraph, '--params', paramsString]
    else:
        cmd = [paths.inflator(), '2D_square', paths.pattern(int(lut.pattern[order[0]]), dim), 
                outPath, '-m',  meshingOptsPath,'--params', paramsString]
    subprocess.check_output(cmd)


#!/usr/bin/env python
import json, sys, paths
import numpy as np
from optimization import *
from LookupTable import LUT, extract as extractLUT
import shutil
import glob
import pdb
import time

if (len(sys.argv) != 4):
    print "usage: ./experiment_targets.py template \"list of Young's module targets\" \"list of patterns\"" 
    print "example: ./horizontal_cover.py template.json \"5 10 20 30\" \"-2 -1 0\""
    sys.exit(-1)

configPath = sys.argv[1]
config = json.load(file(configPath))

EtargetsString = sys.argv[2]
patternsString = sys.argv[3]

Etargets = [float(e) for e in EtargetsString.split(' ')]
patterns = [int(p) for p in patternsString.split(' ')]

print(Etargets)
print(patterns)

folderName = time.strftime("%Y%m%d-%H%M%S")
os.system("mkdir " + folderName)
os.system("cp *.py " + folderName)
os.system("cp *.sh " + folderName)
os.chdir(folderName)

for e in Etargets:
    for p in patterns:
        config['pattern'] = p
        config['targetE'] = e

        tempPath = "temp.txt"
        with open(tempPath, 'w') as outfile:
            json.dump(config, outfile)

        try:
            tableName = "coverage_" + str(p) + "_" + str(e) + ".txt"
            os.system("./horizontal_cover.py " + tempPath)
            #os.system("cp ../horizontalLutTable.txt .")
            os.rename("horizontalLutTable.txt", tableName)
        
            mshCenter = str(p) + "_" + str(e) + "_center" + ".msh"
            mshLeft = str(p) + "_" + str(e) + "_left" + ".msh"
            mshRight = str(p) + "_" + str(e) + "_right" + ".msh"
            os.system("./elasticity_to_mesh.py " + tableName + " " + mshCenter + " " + str(e) + " 0.35")
            os.system("./elasticity_to_mesh.py " + tableName + " " + mshLeft + " " + str(e) + " -1.0")
            os.system("./elasticity_to_mesh.py " + tableName + " " + mshRight + " " + str(e) + " 1.0")

            pngCenter = str(p) + "_" + str(e) + "_center" + ".png"
            pngLeft = str(p) + "_" + str(e) + "_left" + ".png"
            pngRight = str(p) + "_" + str(e) + "_right" + ".png"
            os.system("./msh-to-png.sh " + mshCenter + " " + pngCenter)
            os.system("./msh-to-png.sh " + mshLeft + " " + pngLeft)
            os.system("./msh-to-png.sh " + mshRight + " " + pngRight)

        except:
            print "Run did not work!"

        

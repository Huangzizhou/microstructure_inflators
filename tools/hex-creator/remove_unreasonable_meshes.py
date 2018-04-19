#!/usr/bin/env pythonw
import os
import sys
import fileinput
from shutil import copyfile

if len(sys.argv) < 3:
    print "usage: ./remove_unreasonable_meshes.py <table> <meshes folder>"
    print "example: ./remove_unreasonable_meshes.py results_table1.txt results"
    sys.exit(-1)

dim = 2
volume_fraction = 1.0

# parsing information in Lookup tables and adding them to the chart
i = 0
nu = []
E = []
file_name = []
anisotropy = []

files_to_be_removed = []

tablePath = sys.argv[1]
folder = sys.argv[2]
tableFile = open(tablePath)
for line in tableFile:
    fields = line.strip().split()

    nu.append(float(fields[2]))
    E.append(float(fields[1]))
    anisotropy.append(float(fields[3]))
    file_name.append(fields[5])

    if anisotropy[-1] > 1.05 or anisotropy[-1] < 0.95:
        print "Anisotropy: " + str(anisotropy[-1])
        files_to_be_removed.append(file_name[-1])
    elif (nu[-1] + E[-1]) > 1.0:
        print "Sum: " + str(nu[-1] + E[-1])
        files_to_be_removed.append(file_name[-1])

tableFile.close()

copyfile(tablePath, tablePath + ".backup")

for searched_file in files_to_be_removed:
    for line in fileinput.input(tablePath, inplace=True):
        if not line.find(searched_file) >= 0:
            print line,

for searched_file in files_to_be_removed:
    fixed_name = searched_file[1:]
    os.remove(folder + "/" + fixed_name)

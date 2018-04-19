#!/usr/bin/env python
from __future__ import print_function
import math
import os
import sys
import subprocess

pnorms = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

print("p\tCentered Fin Diff\tDiscrete dJ[v]", end='\n')
for p in pnorms:
    cwd = "../"

    cmd = [cwd + 'ShapeDerivativeValidation_cli', '-b', cwd + "example/dirichlet_thin.json", "-z", cwd + "example/zero_region_empty.json",
            "-m", cwd + "../materials/Russia.material", "--pnorm", str(p), "--usePthRoot", "-o", "validation.tmp", "--perturbationAmplitude",
            "0.0001", "-f", "10.0", cwd + "example/octa_cell_cube_32.msh"]
    #print(cmd)

    app = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    print(str(p), end='\t')

    for line in app.stdout:
        if line.find("Centered difference Stress") != -1:
            lhs, rhs = line.split(":", 1)
            value = float(rhs)
            print(str(value), end='\t')

            #print("Value for centered finite difference is : " + str(value))
        elif line.find("Adjoint discrete shape derivative Stress (delta p)") != -1:
            lhs, rhs = line.split(":", 1)
            value = float(rhs)

            print(str(value), end='\n')
            #print("Value for discrete shape derivative is : " + str(value))


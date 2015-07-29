#!/usr/bin/env 
# Inflates lookup table lines, assuming no parameter constraints.
import os, sys;

if len(sys.argv) != 2:
    print("Usage: inflate.py lut.txt")
    sys.exit(-1)

for i, l in enumerate(file(sys.argv[1])):
    fields = l.split()
    os.system(("{micro}/pattern_optimization/Inflator_cli {micro}/patterns/3D/reference_wires/pattern{pat}.wire"
        + " -I -S0 -p '{params}' -o {mesh}").format(micro=os.environ['MICRO_DIR'],
            pat=fields[0], params=" ".join(fields[5:]), mesh="mesh_%d.msh" % i))

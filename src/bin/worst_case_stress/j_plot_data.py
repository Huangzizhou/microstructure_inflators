#! /usr/bin/env python
# To be run after sample_fields.py, which creates a "samples.hdf" file.
# usage: 
#      ./j_plot_data.py param_name sample_num
import sys
import pandas as pd
import numpy as np

param_name = sys.argv[1]
sample_num = int(sys.argv[2])

data = pd.HDFStore(param_name + '.hdf')
j = data['j']
jdot = data['jdot']
pvalues = data['paramValues']
data.close()

if (sample_num >= j.shape[2]): raise Exception("Invalid sample num");

sname =  j['j'].keys()[sample_num];

jvals = j['j'][sname] 
j_forward = np.diff(jvals) / np.diff(pvalues)
j_centered = np.gradient(jvals, pvalues[1] - pvalues[0])

print "iter\tj\t",
for k in jdot.keys():
    print k,"\t",
print "forward diff\tcentered diff"
for i in range(len(jvals)):
    print "{}\t{}\t".format(pvalues[i], jvals[i]),
    for k in jdot.keys():
        print jdot[k][sname][i],"\t",
    print "{}\t{}".format(j_forward[i] if i < len(j_forward) else 'N/A', j_centered[i])


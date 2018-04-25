#! /usr/bin/env python
# To be run after sample_fields.py, which creates a "samples.hdf" file.
# usage: 
#      ./plot_data.py param_name field_num sample_num
import sys
import pandas as pd
import numpy as np

param_name = sys.argv[1]
field_num = int(sys.argv[2])
sample_num = int(sys.argv[3])

data = pd.HDFStore(param_name + '.hdf')
w = data['w']
wdot = data['wdot']
pvalues = data['paramValues']
data.close()

if (field_num >= w.shape[0]): raise Exception("Invalid field num");
if (sample_num >= w.shape[2]): raise Exception("Invalid sample num");

fname =  w.keys()[field_num];
fdname = wdot.keys()[field_num];
sname =  w[fname].keys()[sample_num];

f = w[fname][sname] 
der = wdot[fdname][sname]
f_forward = np.diff(f) / np.diff(pvalues)
f_centered = np.gradient(f, pvalues[1] - pvalues[0])

print "{}\t{}\t{}\t{}\t{}".format("iter", fname, fdname, "forward diff", "centered diff")
for i in range(len(f)):
    print "{}\t{}\t{}\t{}\t{}".format(pvalues[i], f[i], der[i],
            f_forward[i] if i < len(f_forward) else 'N/A', f_centered[i])

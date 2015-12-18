#!/usr/bin/env python
# usage: lp_hold_wcstress_min_data.py deg area radius
import pandas as pd
import sys

deg,area,radius = sys.argv[1:]

st = pd.HDFStore('d%s_stats.hdf' % deg)

minTable = st['minimum'][area][radius]
argMinTable = st['minimizer'][area][radius]
print "ObjectiveP\tminimizer\tmin"
for k in minTable.keys():
    print "{}\t{}\t{}".format(k, argMinTable[k], minTable[k])

st.close()

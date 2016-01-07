#!/usr/bin/env python
import numpy as np
import pandas as pd
import re, math

for deg in [1,2]:
    minimumData = {}
    minimizerData = {}
    for area in ['0.00004']:
        minimumData[area] = {};
        minimizerData[area] = {};
        for radius in ['0.25', '0.375', '0.5', '0.625', '0.75']:
            minimumData[area][radius] = {}
            minimizerData[area][radius] = {}
            for objectiveP in ['6.0','inf']:
                dataPath = "LpHoleWCStress_12_22_15/d%i_p%s_a%s_r%s.{}.txt" % (deg, objectiveP, area, radius)
                P = 1.0;
                if objectiveP == '6.0':
                    P = 12.0
                    objectiveP = '12.0'

                paramValues, wcObjective = [], []
                for i in range(9):
                    for line in file(dataPath.format(i)):
                        if (re.match("[0-9]+\t[0-9]", line)):
                            fields = line.strip().split()
                            paramValues += [fields[1]]
                            wcObjective += [fields[4]]
                paramValues = np.array(map(float, paramValues))
                wcObjective = np.array(map(float, wcObjective))
                try:
                    minIdx = wcObjective.argmin()
                    minimumData[area][radius][objectiveP] = math.pow(wcObjective[minIdx], 1.0 / P) # Report actual Lp norm, not Lp^p like the optimizer reports
                    minimizerData[area][radius][objectiveP] = paramValues[minIdx]
                except:
                    print "Failed for d%i a%s r%s P%i" % (deg, area, radius, objectiveP)
                    minimumData[area][radius][objectiveP] = float('NaN')
                    minimizerData[area][radius][objectiveP] = float('NaN')

    minPanel    = pd.Panel(minimumData)
    argMinPanel = pd.Panel(minimizerData)

    store = pd.HDFStore('d%i_stats.hdf' % deg)
    store['minimizer'] = argMinPanel
    store['minimum']   = minPanel
    store.close()

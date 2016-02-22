#!/usr/bin/env python
import sys, os, json, math
rootVersion = True
if (len(sys.argv) != 2):
    print "usage: gen_ns_mv_refinement_job.py deg"
    exit()
degree = int(sys.argv[1])
if (degree not in [1, 2]): raise Exception("Degree must be 1 or 2");

WCDir='/home/fjp234/microstructures/worst_case_stress'
expDir='{WCDir}/gradient_validation/lphole'.format(WCDir=WCDir)

cmds = []
for max_vol in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]:
    for nsubdiv in [128,512,2048,8192,16384]:
        for p_center in [1.2, 1.3, 1.4]:
            resultDir='/scratch/fjp234/gradient_validation/wcs/lphole_ns_mvol/d{deg}/ns_{ns}/vol_{vol}/pcenter_{pc}'.format(deg=degree, vol=max_vol, pc=p_center, ns=nsubdiv)
            if not os.path.exists(resultDir): os.makedirs(resultDir)

            for P in range(1, 7):
                for param in range(1):
                    cmds.append({'cmd': 
                        '{WCDir}/GradientComponentValidation {rootFlag}-P{P} --JSWeight 1024 -m /home/fjp234/microstructures/materials/B9Creator.material --inflator lphole --hole_segments {nsubdiv} {expDir}/lphole_job_p{pcenter}.opt -IV -v {max_vol} -d{deg} {param} --range_relative 0.025 10'.format(
                            WCDir=WCDir, expDir=expDir, P=P, param=param, nsubdiv=nsubdiv, max_vol=max_vol, pcenter=p_center, deg=degree,
                            rootFlag=('-R ' if rootVersion else '')),
                    'cwd': resultDir,
                    'stdout': 'P{}_param{}.txt'.format(P, param)})
print json.dumps(cmds, indent=2)

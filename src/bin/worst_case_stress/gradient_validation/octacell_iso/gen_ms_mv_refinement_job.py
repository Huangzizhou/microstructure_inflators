#!/usr/bin/env python
import os, json, sys

if (len(sys.argv) != 2):
    print "usage: gen_ns_mv_refinement_job.py deg"
    exit()
degree = int(sys.argv[1])
if (degree not in [1, 2]): raise Exception("Degree must be 1 or 2");

rootVersion = True
WCDir='/home/fjp234/microstructures/worst_case_stress/'
expDir='{WCDir}/gradient_validation/octacell_iso'.format(WCDir=WCDir)

cmds = []
for max_vol in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]:
    for ms_grid_size in [128,256,512,1024,2048,4096]:
        resultDir='/scratch/fjp234/gradient_validation/wcs/octacell_iso_ms_mvol/d{deg}/ms_{ms}/vol_{vol}'.format(deg=degree, vol=max_vol, ms=ms_grid_size)
        if not os.path.exists(resultDir): os.makedirs(resultDir)
        for P in range(1, 7):
            for param in [0,2,4,5,9]:
                cmds.append({'cmd': 
                    '{WCDir}/GradientComponentValidation {rootFlag}-P{P} --JSWeight 1024 -m /home/fjp234/microstructures/materials/B9Creator.material -p /home/fjp234/microstructures/Luigi/wireinflator2D/meshes/octa_cell.obj {WCDir}/at_target_job_isosurface_octacell.opt -IV -M <({expDir}/mesh_options.sh {mv} {ms}) -d{deg} {param} --range_relative 0.1 10'.format(
                        WCDir=WCDir, expDir=expDir, P=P, param=param, mv=max_vol, ms=ms_grid_size, deg=degree,
                        rootFlag=('-R ' if rootVersion else '')),
                'cwd': resultDir,
                'stdout': 'P{}_param{}.txt'.format(P, param)})
print json.dumps(cmds, indent=2)

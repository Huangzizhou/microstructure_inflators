import json, os
MICRO_DIR=os.environ['MICRO_DIR']
WCDir=MICRO_DIR + '/worst_case_stress'
cmds = []

internal_res = ["1e-4", "5e-5", "1e-5"]
bdry_res = [1024, 2048, 4096, 8192, 16384]
PValues = [1, 1.1, 2, 3, 4, 5, 6, 7, 8]
innerSmoothings = [0.009638, 0.01, 0.02]
params = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13]

for ismooth in innerSmoothings:
    for ir in internal_res:
        for br in bdry_res:
            resultDir=os.environ['SCRATCH'] + '/gradient_validation/wcs/octacell_iso/ismooth{ismooth}/ires{ir}/bres{br}'.format(ismooth=ismooth, ir=ir, br=br)
            if not os.path.exists(resultDir): os.makedirs(resultDir)
            for P in PValues:
                for param in params:
                    cmds.append({
                        'cmd': '{WCDir}/WCSOptimization_cli -R -P{P} --JSWeight 1024 --WCSWeight 1.0  -m {MICRO_DIR}/materials/B9Creator.material -p {MICRO_DIR}/Luigi/wireinflator2D/meshes/octa_cell.obj {WCDir}/gradient_validation_new/octacell/inner_{ismooth}.opt -IV -M <({WCDir}/gradient_validation/octacell_iso/mesh_options.sh {ires} {bres}) -d2 --validateGradientComponent {param} --rangeRelative 0.05 --nsamples 17'.format(
                            WCDir=WCDir, MICRO_DIR=MICRO_DIR, param=param, P=P, ires=ir, bres=br, ismooth=ismooth),
                        'cwd': resultDir,
                        'stdout': 'P{}_param{}.txt'.format(P, param),
                        'stderr': 'P{}_param{}.err.txt'.format(P, param)
                    })
print json.dumps(cmds, indent=2)

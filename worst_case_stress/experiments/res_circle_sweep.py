import json, os
cmd = '$MICRO_DIR/worst_case_stress/WCSRemeshingOptimization_cli -p $MICRO_DIR/worst_case_stress/circle_{subdiv}.msh -m $MICRO_DIR/materials/B9Creator.material $MICRO_DIR/worst_case_stress/at_target_job.opt --JVolWeight {jvolWeight} -d1 -P{P} -n 3000 --maxNormStep 1e-4 -v {maxvol} -O100 -o it_{name} -R'
workingDir = os.environ['SCRATCH'] + '/wcs_bdry_perturb/res_circle_sweep/{subdir}/P{P}/{subdiv}/'
stdout = '{name}.txt'

jobs = []
for P in range(3,24+1):
    for subdiv, mvol in zip([64,90,128,181,256],
                            [0.001, 0.0005, 0.00025, 0.000125, 0.00006125]):
        for weight in [0, 32, 64, 128, 256]:
            jvSweepParameters = dict(
                    subdir="jvol_sweep",
                    name="P{P}_jvol{weight}_sub{sub}".format(P=P, weight=weight, sub=subdiv),
                    P=P,
                    jsWeight=0,
                    maxvol=mvol,
                    jvolWeight=weight,
                    subdiv=subdiv)
            jvWDir = workingDir.format(**jvSweepParameters)
            if (not os.path.exists(jvWDir)): os.makedirs(jvWDir)

            jobs.append({'cmd': cmd.format(**jvSweepParameters), "cwd": workingDir.format(**jvSweepParameters), "stdout": stdout.format(**jvSweepParameters)})

print json.dumps(jobs)

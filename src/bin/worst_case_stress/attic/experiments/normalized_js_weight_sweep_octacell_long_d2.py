import json, os
cmd = '$MICRO_DIR/worst_case_stress/WCSRemeshingOptimization_cli -p $MICRO_DIR/worst_case_stress/octacell_ex_2.msh -m $MICRO_DIR/materials/B9Creator.material $MICRO_DIR/worst_case_stress/at_target_job_octacell_ex_2.opt --JSWeight {jsWeight} --JVolWeight {jvolWeight} -d2 -P{P} -n 6000 --maxNormStep 1e-4 -v 0.00001 -O10 -o it_{name} -R'
workingDir = os.environ['SCRATCH'] + '/wcs_bdry_perturb/normalized_js_weight_sweep_octacell_long_d2/{subdir}'
stdout = '{name}.txt'

jobs = []
for P in [3,5,8,12]:
    for n in range(16):
        weight = 2 ** n
        jsSweepParameters = dict(
                subdir="js_sweep",
                P=P,
                name="P{P}_js{weight}".format(P=P, weight=weight),
                jsWeight=weight,
                jvolWeight=0)
        jsWDir = workingDir.format(**jsSweepParameters)
        if (not os.path.exists(jsWDir)): os.makedirs(jsWDir)

        jobs.append({'cmd': cmd.format(**jsSweepParameters), "cwd": workingDir.format(**jsSweepParameters), "stdout": stdout.format(**jsSweepParameters)})

print json.dumps(jobs)
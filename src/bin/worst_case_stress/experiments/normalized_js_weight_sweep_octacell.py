import json, os
cmd = '$MICRO_DIR/worst_case_stress/WCSRemeshingOptimization_cli -p $MICRO_DIR/worst_case_stress/octacell_ex_2.msh -m $MICRO_DIR/materials/B9Creator.material $MICRO_DIR/worst_case_stress/at_target_job_octacell_ex_2.opt --JSWeight {jsWeight} --JVolWeight {jvolWeight} -d1 -P{P} -n 1500 --maxNormStep 1e-4 -v 0.00001 -O10 -o it_{name} -R'
workingDir = os.environ['SCRATCH'] + '/wcs_bdry_perturb/normalized_js_weight_sweep_octacell/{subdir}'
stdout = '{name}.txt'

jobs = []
for P in [3,5,8,12]:
    for n in range(32):
        weight = 2 ** (n - 16)
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

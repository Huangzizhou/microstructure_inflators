import json, os
cmd = '$MICRO_DIR/worst_case_stress/WCSRemeshingOptimization_cli -p $MICRO_DIR/worst_case_stress/r0.75_p2.msh -m $MICRO_DIR/materials/B9Creator.material $MICRO_DIR/worst_case_stress/at_target_job.opt --JSWeight {jsWeight} --JVolWeight {jvolWeight} -d1 -P{P} -n 5000 --maxNormStep 1e-4 -v 0.001 -O10 -o it_{name}'
workingDir = os.environ['SCRATCH'] + '/wcs_bdry_perturb/{subdir}'
stdout = '{name}.txt'
stderr = '{name}.stderr.txt'

jobs = []
for P in range(3,6+1):
    for n in range(31):
        weight = 2 ** n - 1
        jsSweepParameters = dict(
                subdir="js_sweep",
                P=P,
                name="P{P}_js{weight}".format(P=P, weight=weight),
                jsWeight=weight,
                jvolWeight=0)
        jvSweepParameters = dict(
                subdir="jvol_sweep",
                name="P{P}_jvol{weight}".format(P=P, weight=weight),
                P=P,
                jsWeight=0,
                jvolWeight=weight)
        jobs.append({'cmd': cmd.format(**jsSweepParameters), "cwd": workingDir.format(**jsSweepParameters), "stdout": stdout.format(**jsSweepParameters), "stderr": stderr.format(**jsSweepParameters)})
        jobs.append({'cmd': cmd.format(**jvSweepParameters), "cwd": workingDir.format(**jvSweepParameters), "stdout": stdout.format(**jvSweepParameters), "stderr": stderr.format(**jvSweepParameters)})

print json.dumps(jobs)

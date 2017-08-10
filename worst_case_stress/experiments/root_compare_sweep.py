import json, os, numpy as np
cmd = '$MICRO_DIR/worst_case_stress/WCSRemeshingOptimization_cli -p $MICRO_DIR/worst_case_stress/r0.75_p2.msh -m $MICRO_DIR/materials/B9Creator.material $MICRO_DIR/worst_case_stress/at_target_job.opt -d1 -P{P} -n 5000 --maxNormStep 1e-4 -v 0.001'
workingDir = os.environ['SCRATCH'] + '/wcs_bdry_perturb/root_test/{subdir}'

jobs = []
for P in np.linspace(3, 8, 11):
    baseCmd = cmd.format(P=P)
    stdoutFile = "P{P}.txt".format(P=P)
    stderrFile = "P{P}.stderr.txt".format(P=P)
    jobs.append({'cmd': baseCmd, "cwd": workingDir.format(subdir='no_root'), "stdout": stdoutFile, "stderr": stderrFile})
    jobs.append({'cmd': baseCmd + " -R", "cwd": workingDir.format(subdir='root'), "stdout": stdoutFile, "stderr": stderrFile})

print json.dumps(jobs)

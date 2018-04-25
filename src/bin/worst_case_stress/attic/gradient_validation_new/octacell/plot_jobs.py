import json, os
MICRO_DIR=os.environ['MICRO_DIR']
GVDir=MICRO_DIR + '/worst_case_stress/gradient_validation_new/octacell'
cmds = []

PValues = [1, 1.1, 2, 3, 4, 5, 6, 7, 8]
params = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13]

for P in PValues:
    for param in params:
        cmds.append({
            'cmd': 'python plot.py {} {}'.format(P, param),
            'cwd': GVDir
        })
print json.dumps(cmds, indent=2)

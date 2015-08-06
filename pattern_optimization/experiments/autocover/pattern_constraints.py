import os
def lookup(pat, dim = 3):
    constraints = []
    if (dim == 3):
        for l in file(os.environ['MICRO_DIR'] + "/pattern_optimization/experiments/equality_constraints.txt"):
           fields = l.strip().split('\t')
           if (pat == int(fields[0])):
               constraints = fields[1:]
    return constraints

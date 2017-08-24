#!/usr/bin/env python
import os
import sys
import subprocess
import uuid

if len(sys.argv) != 5:
    print "usage: ./run-optimization.py <wire file> <params file> <volumetric fraction (limit to be considered)> <meshing opts file>"
    print "example: ./run-optimization.py instance.wire instance.params 1.0 refined-meshing_opts.json"
    sys.exit(-1)

wire_path = sys.argv[1]
params_path = sys.argv[2]
vol_frac = float(sys.argv[3])
meshing_file = sys.argv[4]
#target_E = sys.argv[3]
#target_Nu = sys.argv[4]

cwd = os.getcwd()

with open(params_path, 'r') as params_file:
    # read parameters file
    content = params_file.readlines()
    parameters = []
    for line_index in range(0, len(content)):
        for element in content[line_index].split():
            # print element
            parameters.append(float(element))

parameters_string = ' '.join(str(param) for param in parameters)

unique_identifier = str(uuid.uuid1())
folder_path = 'temp-' + unique_identifier
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# Run inflation
out_msh = folder_path + "/out.msh"
cmd = [cwd + '/../../isosurface_inflator/isosurface_cli', '2D_doubly_periodic', wire_path, out_msh, '--params', parameters_string,
       '-m', meshing_file, '-D', 'inflated.msh', '-R', 'replicated.msh']
subprocess.call(cmd)

# Compute elastic properties
table_path = 'temp-' + unique_identifier + '.txt'
cmd = [cwd + '/run-homogenization.py', folder_path, table_path, "../../materials/Russia.material"]
subprocess.call(cmd)

table_file = open(table_path)
for line in table_file:
    fields = line.strip().split()
    experiment = fields[5]

    if experiment == ("#out.msh"):
        initial_E = float(fields[1])
        initial_Nu = float(fields[2])
        break

# Compute targets
E_base = 1.0
nu_base = 0.0
k_base = E_base / (2 * (1 - nu_base))
mu_base = E_base / (2 * (1 + nu_base))

k_top = k_base + (1 - vol_frac) / (vol_frac / (k_base + mu_base) - 1 / k_base)
mu_top = mu_base + (1 - vol_frac) / (vol_frac * (k_base + 2 * mu_base) / (2 * mu_base * (k_base + mu_base)) - 1 / mu_base)
nu_top = (k_top - mu_top) / (k_top + mu_top)
E_top = 2 * mu_top * (1 + nu_top)

target_E = initial_E
target_Nu = (-target_E + E_top) / E_top

#target_Nu = (1 + (initial_Nu - initial_E)) / 2.0
#target_E = -target_Nu + 1.0

print "Original properties: "
print "  Poisson: " + str(initial_Nu)
print "  Young's: " + str(initial_E)

print "Target properties: "
print "  Poisson: " + str(target_Nu)
print "  Young's: " + str(target_E)


# Create job
radiusBounds = [0.00001, 1.0]
offsetBounds = [-0.5, 0.5]
blendingBounds = [0.00001, 0.1]

parameters_string = ', '.join(str(param) for param in parameters)

cmd = ['../../pattern_optimization/GenIsosurfaceJob', wire_path, '-r', str(radiusBounds[0]) + ',' + str(radiusBounds[1]),
    '-o', str(offsetBounds[0]) + ',' + str(offsetBounds[1]), '-b', str(blendingBounds[0]) + ',' + str(blendingBounds[1]),
    '-e', str(target_E) + ',' + str(target_Nu), '--symmetry', '2D_doubly_periodic', '--initialParams', parameters_string]

job_path = os.path.splitext(wire_path)[0] + '.job'
with open(job_path, 'w') as out_job:
    ret = subprocess.call(cmd, stdout=out_job)

# Run optimization
cmd = ['../../worst_case_stress/WCSOptimization_cli', '-p', wire_path, '-m', '../../materials/Russia.material',
       job_path, '-M', meshing_file, '--symmetry', 'doubly_periodic', '--vertexThickness',
       '--WCSWeight', str(0),
       '--JSWeight', str(100.0),
       '--TensorFitConstraint',
       '--proximityRegularizationWeight', str(1.0),
       #'--JIsoWeight', str(1.0),
       '--solver', 'slsqp', '-o', folder_path + '/it', '--deformedCell', '1 0.5 0 0.8660']

out_path = os.path.splitext(wire_path)[0] + '.log'
with open(out_path, 'w') as out_log:
    ret = subprocess.call(cmd, stdout=out_log)

# adding results to table, running homogenization on 'it' instances
for filename in os.listdir(folder_path):
    if filename.startswith("it_"):
        os.rename(folder_path + "/" + filename, folder_path + "/" + filename + ".msh")

table_path = 'temp-' + unique_identifier + '.txt'
cmd = [cwd + '/run-homogenization.py', folder_path, table_path, "../../materials/Russia.material"]
subprocess.call(cmd)

# plotting results
cmd = [cwd + '/plot_luTable.py', 'all', "../../materials/Russia.material", table_path]
subprocess.call(cmd)


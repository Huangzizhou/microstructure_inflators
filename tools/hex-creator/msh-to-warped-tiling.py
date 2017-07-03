#!/usr/bin/env python
import sys, os
sys.path.insert(0, "/Users/davitozoni/Desktop/NYU/Research/Apps/libigl-new/libigl/python")
import pyigl as igl
import numpy as np
from subprocess import call

if len(sys.argv) != 3:
    print "usage: ./msh-to-warped-tiling.py <input mesh file> <output mesh file>"
    print "example: ./msh-to-warped-tiling.py input.msh output.msh"
    sys.exit(-1)

in_path = sys.argv[1]
out_path = sys.argv[2]

V1 = igl.eigen.MatrixXd()
F1 = igl.eigen.MatrixXi()

cwd = os.getcwd()
in_off_path = in_path + '.off'
cmd = [cwd + '/../../../MeshFEM/mesh_convert', in_path, in_off_path]
call(cmd)

igl.readOFF(in_off_path, V1, F1)

npV1 = np.array(V1)
npF1 = np.array(F1)

# translate original microstructure up and down.
# 1 - To the right
npVT = npV1 + [2.0, 0.0, 0.0]
npFT = npF1 + len(npV1)
npV12 = np.concatenate((npV1, npVT), axis=0)
npF12 = np.concatenate((npF1, npFT), axis=0)

# 2 - below
npVT = npV12 + [0.0, 2.0, 0.0]
npFT = npF12 + len(npV12)
npV22 = np.concatenate((npV12, npVT), axis=0)
npF22 = np.concatenate((npF12, npFT), axis=0)

# return to original format (parallelogram)
transformation_matrix = np.array([[1, 0.5, 0.0], [0, 0.8660, 0.0], [0.0, 0.0, 0.0]])
npVC = transformation_matrix.dot(npV22.transpose()).transpose()

VC = igl.eigen.MatrixXd(npVC)
FC = igl.eigen.MatrixXi(npF22)

#viewer = igl.viewer.Viewer()
#viewer.data.clear()
#viewer.data.set_mesh(VC, FC)
#viewer.core.align_camera_center(VC, FC)
#viewer.launch()

out_obj_path = cwd + '/' + out_path + '.obj'
igl.writeOBJ(out_obj_path, VC, FC)
cmd = [cwd + '/../../../MeshFEM/mesh_convert', out_obj_path, cwd + '/' + out_path]
call(cmd)




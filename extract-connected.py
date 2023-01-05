import igl
import numpy as np

v, _, _, f, _, _ = igl.read_obj("debug.obj")

C = igl.connected_components(igl.adjacency_matrix(f))

corner = np.argmin(v, axis=0)[0]
cube_flag = C[1][corner]

if C[0] != 2:
    print("Wrong number of connected components!")
    exit(0)

v_mask = (C[1] != cube_flag)
f_mask = v_mask[f[:,0]]
# for face in f:
#     if v_mask[face[0]]:
#         f_mask.append(True)
#     else:
#         f_mask.append(False)

sv, sf, _, _ = igl.remove_unreferenced(v, f[f_mask,:])
sf[:,[0,1]] = sf[:,[1,0]]
igl.write_obj("debug-clean.obj", sv, sf)
import unittest
import numpy as np

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from MaterialFitterFactory import MaterialFitterFactory
import ResultOutputUtils

class MaterialFitterTest(unittest.TestCase):
    def save_mesh_fields(self, mesh_file, mesh, pressures=None, displacements=None,
            stress_traces=None):
        ResultOutputUtils.save_mesh_fields(mesh_file, mesh, pressures,
                displacements, stress_traces);

    def assertEqualMaterial(self, dim, mat1, mat2):
        ori = np.zeros(dim);
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        self.assertAlmostEqual(
                                mat1.get_material_tensor(i,j,k,l,ori),
                                mat2.get_material_tensor(i,j,k,l,ori));


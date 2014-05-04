import unittest
import numpy as np

import LinearElasticitySettings
from Material import Material
from mesh_io import load_mesh, save_mesh
from MaterialFitterFactory import MaterialFitterFactory
from MaterialFitterTest import MaterialFitterTest

import PyAssembler

class MaterialFitter3DTest(MaterialFitterTest):
    def setUp(self):
        self.material_file = "examples/symmetric_3D.material";
        self.material = Material(3, self.material_file).material;
        self.mesh= load_mesh("examples/solid_cube.msh");
        factory = MaterialFitterFactory(self.mesh, self.material_file);
        self.fitter = factory.create("symmetric");

    def test_3D(self):
        fitter = self.fitter;
        fitter.fit();
        self.assertAlmostEqual(0.0, fitter.residual_error);

        homogenized_material = PyAssembler.Material.create_symmetric(
                self.material.get_density(), fitter.elasticity_tensor);
        self.assertEqualMaterial(3, self.material, homogenized_material);


import unittest
import numpy as np

import LinearElasticitySettings
from Material import Material
from mesh_io import load_mesh, save_mesh
from MaterialFitterFactory import MaterialFitterFactory
from MaterialFitterTest import MaterialFitterTest

import PyAssembler

class OrthotropicMaterialFitter3DTest(MaterialFitterTest):
    def setUp(self):
        self.mesh= load_mesh("examples/solid_cube.msh");
        factory = MaterialFitterFactory(self.mesh, None);
        self.fitter = factory.create("orthotropic");

    def test_3D(self):
        fitter = self.fitter;
        fitter.fit();
        self.save_mesh_fields("tmp3D.msh", fitter.mesh,
                fitter.pressures,
                fitter.displacements,
                fitter.stress_traces);
        self.save_mesh_fields("tmp3D_coarse.msh", fitter.coarse_mesh,
                displacements = fitter.coarse_displacements);
        self.assertEqual(3, len(fitter.youngs_modulus));
        self.assertEqual(6, len(fitter.poisson_ratio));
        self.assertEqual(3, len(fitter.shear_modulus));
        self.assertAlmostEqual(1.0, np.amax(fitter.youngs_modulus));
        self.assertAlmostEqual(1.0, np.amin(fitter.youngs_modulus));
        self.assertAlmostEqual(0.0, np.amax(fitter.poisson_ratio));
        self.assertAlmostEqual(0.0, np.amin(fitter.poisson_ratio));
        self.assertAlmostEqual(0.5, np.amax(fitter.shear_modulus));
        self.assertAlmostEqual(0.5, np.amin(fitter.shear_modulus));

    def test_orthotropic_cube(self):
        self.mesh= load_mesh("examples/solid_cube.msh");
        factory = MaterialFitterFactory(self.mesh,
                "examples/orthotropic.material");
        fitter = factory.create("orthotropic");
        fitter.fit();

        material = Material(3, "examples/orthotropic.material").material;
        homogenized_material = PyAssembler.Material.create_symmetric(
                material.get_density(), fitter.elasticity_tensor);

        self.assertEqualMaterial(3, material, homogenized_material);


import unittest
import numpy as np

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from MaterialFitterFactory import MaterialFitterFactory
from MaterialFitterTest import MaterialFitterTest

class OrthotropicMaterialFitter2DTest(MaterialFitterTest):
    def setUp(self):
        self.mesh= load_mesh("examples/solid_square.obj");
        factory = MaterialFitterFactory(self.mesh,
                "examples/isotropic.material");
        self.fitter = factory.create("orthotropic");

    def test_2D(self):
        fitter = self.fitter;
        fitter.fit();
        self.save_mesh_fields("tmp2D.msh", fitter.mesh,
                fitter.pressures,
                fitter.displacements,
                fitter.stress_traces);
        self.save_mesh_fields("tmp2D_coarse.msh", fitter.coarse_mesh,
                displacements = fitter.coarse_displacements);
        self.assertEqual(2, len(fitter.youngs_modulus));
        self.assertEqual(2, len(fitter.poisson_ratio));
        self.assertEqual(1, len(fitter.shear_modulus));
        self.assertAlmostEqual(1.0, fitter.youngs_modulus[0]);
        self.assertAlmostEqual(1.0, fitter.youngs_modulus[1]);
        self.assertAlmostEqual(0.0, fitter.poisson_ratio[0]);
        self.assertAlmostEqual(0.0, fitter.poisson_ratio[1]);
        self.assertAlmostEqual(0.5, fitter.shear_modulus[0]);


import unittest

import LinearElasticitySettings
from mesh_io import load_mesh
from IsotropicMaterialOptimizer import IsotropicMaterialOptimizer
from timethis import timethis

class IsotropicMaterialOptimizerTest(unittest.TestCase):
    def setUp(self):
        self.mesh = load_mesh("example/square_8x8.obj");
        self.optimizer = IsotropicMaterialOptimizer(self.mesh);

    def test_initialization(self):
        self.assertTrue(self.mesh.has_attribute(self.optimizer.young_field_name));
        self.assertTrue(self.mesh.has_attribute(self.optimizer.poisson_field_name));
        self.assertTrue(self.mesh.has_attribute(self.optimizer.grad_young_field_name));
        self.assertTrue(self.mesh.has_attribute(self.optimizer.grad_poisson_field_name));
        timethis.summarize();


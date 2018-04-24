import numpy as np
import unittest

import LinearElasticitySettings
from mesh_io import save_mesh, load_mesh

from BoxMeshGenerator import generate_box_mesh

class BoxMeshGeneratorTest(unittest.TestCase):
    def test_2D_mesh(self):
        box_min = np.zeros(2);
        box_max = np.ones(2);
        num_samples = 2;
        mesh, quad_indices = generate_box_mesh(box_min, box_max, num_samples);
        self.assertEqual((num_samples+1)**2, mesh.num_vertices);
        self.assertEqual(num_samples**2*2, mesh.num_faces);

    def test_3D_mesh(self):
        box_min = np.zeros(3);
        box_max = np.ones(3)*10;
        num_samples = 2;
        mesh, hex_indices = generate_box_mesh(box_min, box_max, num_samples);
        #self.assertEqual((num_samples+1)**3, mesh.num_vertices);
        #self.assertEqual(num_samples**2*2*6, mesh.num_faces);
        save_mesh("box.msh", mesh);


import unittest
from Triangulation import Triangulation

import numpy as np
from numpy.linalg import norm, det

class TriangulationTest(unittest.TestCase):
    def test_triangle(self):
        vertices = np.eye(3);
        faces = np.array([0, 1, 2], dtype=int).reshape((1, 3));
        project_dir = np.ones(3);

        tri = Triangulation(vertices, faces, project_dir);
        tri.triangulate(0.1);

        self.assertAlmostEqual(1.0, det(tri.coordinate_frame));
        intercept = np.dot(vertices[0], project_dir);
        for v in tri.vertices:
            self.assertAlmostEqual(intercept, np.dot(v, project_dir));

    def test_quad(self):
        vertices = np.array([
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0] ]);
        faces = np.array([
            [0, 1, 2],
            [0, 2, 3] ]);
        project_dir = np.array([0.0, 0.0, 1.0]);

        tri = Triangulation(vertices, faces, project_dir);
        tri.triangulate(1.0);

        self.assertEqual(2, len(tri.faces));
        self.assertEqual(4, len(tri.vertices));

        self.assertAlmostEqual(1.0, det(tri.coordinate_frame));
        intercept = np.dot(vertices[0], project_dir);
        for v in tri.vertices:
            self.assertAlmostEqual(intercept, np.dot(v, project_dir));


import unittest
from Triangulation import Triangulation

import numpy as np
from numpy.linalg import norm, det

class TriangulationTest(unittest.TestCase):
    def assert_triangulation_is_valid(self, tri, proj_dir):
        vertices = tri.vertices;
        faces = tri.faces;
        for face in faces:
            self.assertEqual(3, len(face));
            v0 = vertices[face[0]];
            v1 = vertices[face[1]];
            v2 = vertices[face[2]];
            n = np.cross(v1 - v0, v2 - v0);
            n /= norm(n);
            proj_dir /= norm(proj_dir);
            self.assertLess(0.5, n.dot(proj_dir));

    def test_triangle(self):
        vertices = np.eye(3);
        faces = np.array([0, 1, 2], dtype=int).reshape((1, 3));
        project_dir = np.ones(3);

        tri = Triangulation(vertices, faces, project_dir);
        tri.triangulate(0.1);

        intercept = np.dot(vertices[0], project_dir);
        for v in tri.vertices:
            self.assertAlmostEqual(intercept, np.dot(v, project_dir));
        self.assert_triangulation_is_valid(tri, project_dir);

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

        intercept = np.dot(vertices[0], project_dir);
        for v in tri.vertices:
            self.assertAlmostEqual(intercept, np.dot(v, project_dir));
        self.assert_triangulation_is_valid(tri, project_dir);


import unittest
from BoxIntersection import BoxIntersection

import numpy as np
from numpy.linalg import norm

class BoxIntersectionTest(unittest.TestCase):
    def setUp(self):
        bbox_min = np.ones(3) * -1;
        bbox_max = np.ones(3);
        self.box = BoxIntersection(bbox_min, bbox_max);

    def distance(self, p, points):
        dist = norm(points - p, axis=1);
        return np.amin(dist);

    def assert_coplanar(self, p, p1, p2, p3):
        normal = np.cross(p2 - p1, p3 - p1);
        height = np.dot(normal, p - p1);
        self.assertAlmostEqual(0.0, height);

    def assert_inside(self, p, p1, p2, p3):
        eps = 1e-6;
        normal = np.cross(p2 - p1, p3 - p1);
        e1 = p1 - p;
        e2 = p2 - p;
        e3 = p3 - p;
        a1 = np.dot(np.cross(e1, e2), normal);
        a2 = np.dot(np.cross(e2, e3), normal);
        a3 = np.dot(np.cross(e3, e1), normal);
        area = np.array([a1, a2, a3]);
        self.assertTrue(np.all(area > -eps));

    def assert_intersection_is_valid(self, intersection, p1, p2, p3):
        for p in intersection:
            self.assert_inside(p, p1, p2, p3);
            self.assert_coplanar(p, p1, p2, p3);

    def test_not_intersect(self):
        p0 = np.zeros(3);
        p1 = np.array([0.1, 0.0, 0.0]);
        p2 = np.array([0.0, 0.1, 0.0]);
        vertices = np.vstack([p0, p1, p2]);
        faces = np.array([[0, 1, 2]]);
        intersection = self.box.intersect(vertices, faces);
        self.assertEqual(0, len(intersection));

    def test_single_face_intersect(self):
        p0 = np.zeros(3);
        p1 = np.array([2.0, 0.0, 0.0]);
        p2 = np.array([0.0, 0.5, 0.0]);
        vertices = np.vstack([p0, p1, p2]);
        faces = np.array([[0, 1, 2]]);
        intersection = self.box.intersect(vertices, faces);
        self.assertEqual(2, len(intersection));
        self.assert_intersection_is_valid(intersection, p0, p1, p2);

    def test_double_face_intersect(self):
        p0 = np.zeros(3);
        p1 = np.array([1.6, 0.0, 0.0]);
        p2 = np.array([0.0, 1.6, 0.0]);
        vertices = np.vstack([p0, p1, p2]);
        faces = np.array([[0, 1, 2]]);
        intersection = self.box.intersect(vertices, faces);
        self.assertEqual(4, len(intersection));
        self.assert_intersection_is_valid(intersection, p0, p1, p2);

        expected_intersections = np.array([
                [1.0, 0.0, 0.0],
                [1.0, 0.6, 0.0],
                [0.0, 1.0, 0.0],
                [0.6, 1.0, 0.0] ]);

        for p in expected_intersections:
            self.assertAlmostEqual(0.0, self.distance(p, intersection));

    def test_double_face_intersect_2(self):
        p0 = np.zeros(3);
        p1 = np.array([2.6, 0.0, 0.0]);
        p2 = np.array([0.0, 2.6, 0.0]);
        vertices = np.vstack([p0, p1, p2]);
        faces = np.array([[0, 1, 2]]);
        intersection = self.box.intersect(vertices, faces);
        self.assertEqual(3, len(intersection));
        self.assert_intersection_is_valid(intersection, p0, p1, p2);

        expected_intersections = np.array([
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0] ]);

        for p in expected_intersections:
            self.assertAlmostEqual(0.0, self.distance(p, intersection));

    def test_triple_face_intersect_2(self):
        p0 = np.array([0.4, 1.1, 1.1]);
        p1 = np.array([1.1, 0.4, 1.1]);
        p2 = np.array([1.1, 1.1, 0.4]);
        vertices = np.vstack([p0, p1, p2]);
        faces = np.array([[0, 1, 2]]);
        intersection = self.box.intersect(vertices, faces);
        self.assertEqual(3, len(intersection));
        self.assert_intersection_is_valid(intersection, p0, p1, p2);

        expected_intersections = np.array([
                [0.6, 1.0, 1.0],
                [1.0, 0.6, 1.0],
                [1.0, 1.0, 0.6] ]);

        for p in expected_intersections:
            self.assertAlmostEqual(0.0, self.distance(p, intersection));

    def test_uniqueness(self):
        p0 = np.array([0.0, 0.0, 1.0]);
        p1 = np.array([0.2, 0.2, 1.0]);
        p2 = np.array([0.2, 0.0, 1.0]);
        p3 = np.array([0.0, 0.0, 1.001]);
        p4 = np.array([0.2, 0.2, 1.001]);
        p5 = np.array([0.2, 0.0, 1.001]);
        vertices = np.vstack([p0, p1, p2, p3, p4, p5]);
        faces = np.array([
            [0, 1, 2],
            [3, 4, 5],
            [0, 1, 3],
            [3, 1, 4],
            [1, 2, 4],
            [4, 2, 5],
            [2, 0, 5],
            [5, 0, 3] ]);
        intersection = self.box.intersect(vertices, faces);
        self.assertEqual(3, len(intersection));
        


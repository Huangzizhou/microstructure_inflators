import unittest
import numpy as np
from WirePattern import WirePattern
from WireInflator import WireInflator

class WirePatternTest(unittest.TestCase):
    def setUp(self):
        self.vertices = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                ];
        self.edges = [
                [0, 1],
                [0, 2],
                [0, 3] ];

    def test_creation(self):
        pattern = WirePattern();
        pattern.set_single_cell(self.vertices, self.edges);
        bbox_max = np.amax(self.vertices, axis=0);
        bbox_min = np.amin(self.vertices, axis=0);
        centroid = 0.5 * (bbox_max + bbox_min);

        self.assertListEqual((bbox_min - centroid).tolist(),
                pattern.pattern_bbox_min.tolist());
        self.assertListEqual((bbox_max - centroid).tolist(),
                pattern.pattern_bbox_max.tolist());
        self.assertListEqual((bbox_max - bbox_min).tolist(),
                pattern.pattern_bbox_size.tolist());

    def test_tile(self):
        pattern = WirePattern();
        pattern.set_single_cell(self.vertices, self.edges);
        pattern.tile([2, 2, 2]);
        bbox_min = np.amin(pattern.wire_vertices, axis=0);
        bbox_max = np.amax(pattern.wire_vertices, axis=0);
        self.assertListEqual([-1.0, -1.0, -1.0], bbox_min.tolist());
        self.assertListEqual([ 1.0,  1.0,  1.0], bbox_max.tolist());

    def test_remove_duplicates(self):
        pattern = WirePattern();
        pattern.set_single_cell(self.vertices, self.edges);
        pattern.tile([2, 2, 2]);
        self.assertEqual(20, len(pattern.wire_vertices));
        inflator = WireInflator(pattern.wire_network);
        inflator.inflate(0.1);
        inflator.save("tmp2.obj");


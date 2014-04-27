import unittest
from WireNetwork import WireNetwork
from PeriodicWireInflator import PeriodicWireInflator

import numpy as np
from numpy.linalg import norm

class PeriodicWireInflatorTest(unittest.TestCase):
    def setUp(self):
        self.wire_file = "examples/cube.wire";
        self.wire_network = WireNetwork();
        self.wire_network.load_from_file(self.wire_file);

    def test_create(self):
        thickness = 0.1
        inflator = PeriodicWireInflator(self.wire_network);
        inflator.inflate(thickness);
        inflator.save("tmp3.obj");
        bbox_min, bbox_max = self.wire_network.bbox;
        for v in inflator.mesh_vertices:
            self.assertTrue(np.all(v < bbox_max + 4*thickness));
            self.assertTrue(np.all(v > bbox_min - 4*thickness));


    def test_edge_loop(self):
        thickness = 0.1;
        inflator = PeriodicWireInflator(self.wire_network);
        inflator.inflate(thickness);
        for i,edge in enumerate(inflator.wire_network.edges):
            loop_1 = inflator.edge_loops[i, 0, :, :];
            loop_2 = inflator.edge_loops[i, 1, :, :];
            center_1 = np.mean(loop_1, axis=0);
            center_2 = np.mean(loop_2, axis=0);
            dist_1 = norm(center_1 - inflator.wire_network.vertices[edge[0]]);
            dist_2 = norm(center_2 - inflator.wire_network.vertices[edge[1]]);
            self.assertLess(dist_1, thickness*4);
            self.assertLess(dist_2, thickness*4);

    def test_edge_pipe(self):
        eps = 1e-3
        thickness = 0.1;
        inflator = PeriodicWireInflator(self.wire_network);
        inflator.inflate(thickness);
        bbox_min, bbox_max = self.wire_network.bbox;
        for i,edge in enumerate(inflator.wire_network.edges):
            if np.any(inflator.phantom_vertex_map[edge] < 0):
                continue;
            loop_1 = inflator.edge_loops[i, 0, :, :];
            loop_2 = inflator.edge_loops[i, 1, :, :];
            center_1 = np.mean(loop_1, axis=0);
            center_2 = np.mean(loop_2, axis=0);

            self.assertTrue(np.any(center_1 < bbox_max + eps));
            self.assertTrue(np.any(center_2 < bbox_max + eps));
            self.assertTrue(np.any(center_1 > bbox_min - eps));
            self.assertTrue(np.any(center_2 > bbox_min - eps));


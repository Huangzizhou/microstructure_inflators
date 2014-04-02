import unittest
from WireNetwork import WireNetwork
from WireInflator import WireInflator

import numpy as np
from numpy.linalg import norm

class WireInflatorTest(unittest.TestCase):
    def setUp(self):
        self.wire_file = "examples/example.wire";
        self.wire_network = WireNetwork();
        self.wire_network.load_from_file(self.wire_file);

    def test_create(self):
        inflator = WireInflator(self.wire_network);
        inflator.inflate(0.1);
        inflator.save("tmp.obj");
        self.assertEqual(len(self.wire_network.edges), inflator.edge_loops.shape[0]);

    def test_edge_loops(self):
        thickness = 0.1;
        inflator = WireInflator(self.wire_network);
        inflator.inflate(thickness);
        for i,edge in enumerate(self.wire_network.edges):
            loop_1 = inflator.edge_loops[i, 0, :, :];
            loop_2 = inflator.edge_loops[i, 1, :, :];
            center_1 = np.mean(loop_1, axis=0);
            center_2 = np.mean(loop_2, axis=0);
            dist_1 = norm(center_1 - self.wire_network.vertices[edge[0]]);
            dist_2 = norm(center_2 - self.wire_network.vertices[edge[1]]);
            self.assertLess(dist_1, thickness*4);
            self.assertLess(dist_2, thickness*4);

    def test_edge_loop_indices(self):
        thickness = 0.1;
        inflator = WireInflator(self.wire_network);
        inflator.inflate(thickness, clean_up=False);
        for i, edge in enumerate(self.wire_network.edges):
            loop_1 = inflator.edge_loops[i, 0, :, :];
            loop_2 = inflator.edge_loops[i, 1, :, :];
            loop_1_indices = inflator.edge_loop_indices[i*2];
            loop_2_indices = inflator.edge_loop_indices[i*2+1];
            for i in range(4):
                v1 = inflator.mesh_vertices[loop_1_indices[i], :]
                v2 = inflator.mesh_vertices[loop_2_indices[i], :]
                self.assertListEqual(v1.tolist(), loop_1[i, :].tolist());
                self.assertListEqual(v2.tolist(), loop_2[i, :].tolist());

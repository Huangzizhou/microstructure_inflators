import unittest
from WireNetwork import WireNetwork
from math import sqrt

class WireNetworkTest(unittest.TestCase):
    def setUp(self):
        self.wire_file = "examples/example.wire"
        self.wire_frame = WireNetwork();
        self.wire_frame.load_from_file(self.wire_file);

    def test_creation(self):
        self.assertEqual(3, self.wire_frame.num_vertices);
        self.assertEqual(3, self.wire_frame.num_edges);

    def test_connectivity(self):
        self.assertEqual([1,2], sorted(self.wire_frame.vertex_neighbors[0]));
        self.assertEqual([0,2], sorted(self.wire_frame.vertex_neighbors[1]));
        self.assertEqual([0,1], sorted(self.wire_frame.vertex_neighbors[2]));

    def test_total_length(self):
        self.assertAlmostEqual(2.0+sqrt(2), self.wire_frame.total_wire_length);


import unittest
import numpy as np
from WireNetwork import WireNetwork
from utils.find_file import find_file
from math import sqrt

class WireNetworkTest(unittest.TestCase):
    def setUp(self):
        self.load_wire("examples/example.wire");

    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        self.wire_frame = WireNetwork();
        self.wire_frame.load_from_file(wire_file);

    def test_creation(self):
        self.assertEqual(3, self.wire_frame.num_vertices);
        self.assertEqual(3, self.wire_frame.num_edges);

    def test_connectivity(self):
        self.assertEqual([1,2], sorted(self.wire_frame.vertex_neighbors[0]));
        self.assertEqual([0,2], sorted(self.wire_frame.vertex_neighbors[1]));
        self.assertEqual([0,1], sorted(self.wire_frame.vertex_neighbors[2]));

    def test_total_length(self):
        self.assertAlmostEqual(2.0+sqrt(2), self.wire_frame.total_wire_length);

    def test_trim(self):
        self.load_wire("examples/example2.wire");
        self.assertEqual(20, self.wire_frame.num_vertices);
        self.assertEqual(30, self.wire_frame.num_edges);

        self.wire_frame.trim();

        self.assertEqual(14, self.wire_frame.num_vertices);
        self.assertEqual(24, self.wire_frame.num_edges);

    def test_2D_wire(self):
        self.load_wire("examples/patterns/box_2D.wire");
        self.assertEqual(2, self.wire_frame.dim);
        self.assertEqual(12, self.wire_frame.num_vertices);
        self.assertEqual(12, self.wire_frame.num_edges);

        bbox_center = self.wire_frame.bbox_center;
        self.assertListEqual([0.5, 0.5], bbox_center.tolist());

    def test_symmetry_orbits(self):
        self.load_wire("patterns/3D/brick5.wire");
        self.wire_frame.compute_symmetry_orbits();
        vertex_orbits = self.wire_frame.attributes["orthotropic_symmetry_vertex_orbit"];
        edge_orbits = self.wire_frame.attributes["orthotropic_symmetry_edge_orbit"];

        self.assertEqual(len(vertex_orbits), self.wire_frame.num_vertices);
        self.assertEqual(len(edge_orbits), self.wire_frame.num_edges);

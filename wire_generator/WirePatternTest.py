import unittest
import numpy as np
from WireNetwork import WireNetwork
from WirePattern import WirePattern

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

    def test_vertex_attributes(self):
        cell = WireNetwork();
        cell.load_from_file("examples/cube.wire");
        cell.attributes.add("symmetry_orbit");
        cell_orbits = cell.attributes["symmetry_orbit"];

        pattern = WirePattern();
        pattern.set_single_cell_from_wire_network(cell);
        pattern.tile([2, 2, 2]);
        tiled_wires = pattern.wire_network;

        self.assertTrue("symmetry_orbit" in tiled_wires.attributes);

        orbits = tiled_wires.attributes["symmetry_orbit"];
        self.assertEqual(len(tiled_wires.vertices), len(orbits));
        self.assertSetEqual(set(cell_orbits), set(orbits));

        for i,vi in enumerate(pattern.pattern_vertex_map):
            self.assertEqual(orbits[i], cell_orbits[vi]);


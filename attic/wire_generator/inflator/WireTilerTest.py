import unittest
import numpy as np
from core.WireNetwork import WireNetwork
from WireTiler import WireTiler

import PyMesh

class WireTilerTest(unittest.TestCase):
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

    def load_mesh(self, mesh_file):
        factory = PyMesh.MeshFactory();
        factory.load_file(mesh_file);
        mesh = factory.create();
        return mesh;

    def test_creation(self):
        pattern = WireTiler();
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
        pattern = WireTiler();
        pattern.set_single_cell(self.vertices, self.edges);
        pattern.tile([2, 2, 2]);
        bbox_min = np.amin(pattern.wire_vertices, axis=0);
        bbox_max = np.amax(pattern.wire_vertices, axis=0);
        self.assertListEqual([-1.0, -1.0, -1.0], bbox_min.tolist());
        self.assertListEqual([ 1.0,  1.0,  1.0], bbox_max.tolist());

    def test_remove_duplicates(self):
        pattern = WireTiler();
        pattern.set_single_cell(self.vertices, self.edges);
        pattern.tile([2, 2, 2]);
        self.assertEqual(20, len(pattern.wire_vertices));

    def test_vertex_attributes(self):
        cell = WireNetwork();
        cell.load_from_file("examples/cube.wire");
        cell.attributes.add("orthotropic_symmetry_vertex_orbit");
        cell_orbits = cell.attributes["orthotropic_symmetry_vertex_orbit"];

        pattern = WireTiler();
        pattern.set_single_cell_from_wire_network(cell);
        pattern.tile([2, 2, 2]);
        tiled_wires = pattern.wire_network;

        self.assertTrue("orthotropic_symmetry_vertex_orbit" in tiled_wires.attributes);

        orbits = tiled_wires.attributes["orthotropic_symmetry_vertex_orbit"];
        self.assertEqual(len(tiled_wires.vertices), len(orbits));
        self.assertListEqual(
                np.unique(cell_orbits.ravel()).tolist(),
                np.unique(orbits.ravel()).tolist() );

    def test_tile_hex(self):
        cell = WireNetwork();
        cell.load_from_file("examples/cube.wire");
        mesh = self.load_mesh("examples/bar.msh");

        pattern = WireTiler();
        pattern.set_single_cell_from_wire_network(cell);
        pattern.tile_hex_mesh(mesh);

        tiled_wires = pattern.wire_network;

        wire_bbox_min, wire_bbox_max = tiled_wires.bbox;
        wire_bbox_size = wire_bbox_max - wire_bbox_min;
        mesh_bbox_min = np.amin(mesh.get_vertices().reshape((-1, 3)), axis=0);
        mesh_bbox_max = np.amax(mesh.get_vertices().reshape((-1, 3)), axis=0);
        mesh_bbox_size = mesh_bbox_max - mesh_bbox_min;
        self.assertListEqual(mesh_bbox_size.tolist(), wire_bbox_size.tolist());


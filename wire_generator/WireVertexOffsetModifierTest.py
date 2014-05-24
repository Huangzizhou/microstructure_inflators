import unittest
import copy
import numpy as np
from numpy.linalg import norm
import json

from WireVertexOffsetModifier import WireVertexOffsetModifier
from WireNetwork import WireNetwork

class WireVertexOffsetModifierTest(unittest.TestCase):
    def setUp(self):
        self.config = {
                "type": "vertex_orbit",
                "effective_orbits": [0],
                "offset_percentages": [[0.5, 0.5, 0.5]]
                };

    def load_wire(self, filename):
        wire_network = WireNetwork();
        wire_network.load_from_file(filename);
        return wire_network;

    def load_orbit(self, orbit_file):
        with open(orbit_file, 'r') as fin:
            orbit_config = json.load(fin);
        return orbit_config;

    def compute_orbit_centers(self, wire_network, orbits):
        centers = [];
        for orbit in orbits:
            vertices = wire_network.vertices[orbit];
            bbox_min = np.amin(vertices, axis=0);
            bbox_max = np.amax(vertices, axis=0);
            center = 0.5 * (bbox_min + bbox_max);
            centers.append(center);
        return centers;

    def test_creation(self):
        orbit_0_percentages = [0.1, 0.2, 0.3];
        wire_network = self.load_wire("examples/cube.wire");
        ori_wire_network = copy.deepcopy(wire_network);
        self.config["orbit_file"] = "examples/cube.orbit"
        self.config["offset_percentages"] = [orbit_0_percentages];

        orbits = self.load_orbit(self.config["orbit_file"]);
        vertex_orbits = orbits["vertex_orbits"];

        orbit_centers = self.compute_orbit_centers(wire_network,
                vertex_orbits);

        modifier = WireVertexOffsetModifier(self.config);
        modifier.modify(wire_network);

        for i in range(wire_network.num_vertices):
            if i in vertex_orbits[0]:
                orbit_center = orbit_centers[0];
                ori_v = ori_wire_network.vertices[i];
                mod_v = wire_network.vertices[i];
                ratios = (mod_v - ori_v) / (ori_v - orbit_center);
                self.assertAlmostEqual(orbit_0_percentages[0], ratios[0]);
                self.assertAlmostEqual(orbit_0_percentages[1], ratios[1]);
                self.assertAlmostEqual(orbit_0_percentages[2], ratios[2]);
            else:
                self.assertListEqual(
                        ori_wire_network.vertices[i].tolist(),
                        wire_network.vertices[i].tolist());

    def test_formula(self):
        x = 0.1;
        y = 0.5;
        orbit_0_percentages = ["{x}", "{y}", "{x} + {y}"];
        wire_network = self.load_wire("examples/cube.wire");
        ori_wire_network = copy.deepcopy(wire_network);
        self.config["orbit_file"] = "examples/cube.orbit"
        self.config["offset_percentages"] = [orbit_0_percentages];

        orbits = self.load_orbit(self.config["orbit_file"]);
        vertex_orbits = orbits["vertex_orbits"];

        orbit_centers = self.compute_orbit_centers(wire_network,
                vertex_orbits);

        modifier = WireVertexOffsetModifier(self.config);
        modifier.modify(wire_network, x=x, y=y);

        for i in range(wire_network.num_vertices):
            if i in vertex_orbits[0]:
                orbit_center = orbit_centers[0];
                ori_v = ori_wire_network.vertices[i];
                mod_v = wire_network.vertices[i];
                ratios = (mod_v - ori_v) / (ori_v - orbit_center);
                self.assertAlmostEqual(x, ratios[0]);
                self.assertAlmostEqual(y, ratios[1]);
                self.assertAlmostEqual(x+y, ratios[2]);
            else:
                self.assertListEqual(
                        ori_wire_network.vertices[i].tolist(),
                        wire_network.vertices[i].tolist());


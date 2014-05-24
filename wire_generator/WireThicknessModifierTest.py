import unittest
import numpy as np
import json

from WireThicknessModifier import WireThicknessModifier
from WireNetwork import WireNetwork

class WireThicknessModifierTest(unittest.TestCase):
    def setUp(self):
        self.config = {
                "type": "vertex_orbit",
                "effective_orbits": [0],
                "thickness": [1.0],
                "default": 0.5
                };

    def load_wire(self, filename):
        wire_network = WireNetwork();
        wire_network.load_from_file(filename);
        return wire_network;

    def load_orbit(self, orbit_file):
        with open(orbit_file, 'r') as fin:
            orbit_config = json.load(fin);
        return orbit_config;

    def test_creation(self):
        wire_network = self.load_wire("examples/cube.wire");
        self.config["orbit_file"] = "examples/cube.orbit"
        orbits = self.load_orbit(self.config["orbit_file"]);

        modifier = WireThicknessModifier(self.config);
        modifier.modify(wire_network);

        self.assertTrue("vertex_thickness" in wire_network.attributes);
        thicknesses = wire_network.attributes["vertex_thickness"];

        self.assertEqual(wire_network.num_vertices, len(thicknesses));
        for i,thickness in enumerate(thicknesses):
            if i in orbits["vertex_orbits"][0]:
                self.assertAlmostEqual(1.0, thickness);
            else:
                self.assertAlmostEqual(0.5, thickness);

    def test_creation2(self):
        wire_network = self.load_wire("examples/example2.wire");
        self.config["orbit_file"] = "examples/example2.orbit"
        orbits = self.load_orbit(self.config["orbit_file"]);

        modifier = WireThicknessModifier(self.config);
        modifier.modify(wire_network);

        self.assertTrue("vertex_thickness" in wire_network.attributes);
        thicknesses = wire_network.attributes["vertex_thickness"];

        self.assertEqual(wire_network.num_vertices, len(thicknesses));
        for i,thickness in enumerate(thicknesses):
            if i in orbits["vertex_orbits"][0]:
                self.assertAlmostEqual(1.0, thickness);
            else:
                self.assertAlmostEqual(0.5, thickness);

    def test_creation_edge_orbits(self):
        wire_network = self.load_wire("examples/example2.wire");
        self.config["orbit_file"] = "examples/example2.orbit";
        self.config["type"] = "edge_orbit";
        orbits = self.load_orbit(self.config["orbit_file"]);

        modifier = WireThicknessModifier(self.config);
        modifier.modify(wire_network);

        self.assertTrue("edge_thickness" in wire_network.attributes);
        thicknesses = wire_network.attributes["edge_thickness"];

        self.assertEqual(wire_network.num_edges, len(thicknesses));
        for i,thickness in enumerate(thicknesses):
            if i in orbits["edge_orbits"][0]:
                self.assertAlmostEqual(1.0, thickness);
            else:
                self.assertAlmostEqual(0.5, thickness);


    def test_formula(self):
        wire_network = self.load_wire("examples/example2.wire");
        self.config["orbit_file"] = "examples/example2.orbit";
        self.config["thickness"] = ["{x} + {y}"];
        orbits = self.load_orbit(self.config["orbit_file"]);

        modifier = WireThicknessModifier(self.config);
        modifier.modify(wire_network, x=10, y=1);

        self.assertTrue("vertex_thickness" in wire_network.attributes);
        thicknesses = wire_network.attributes["vertex_thickness"];

        self.assertEqual(wire_network.num_vertices, len(thicknesses));
        for i,thickness in enumerate(thicknesses):
            if i in orbits["vertex_orbits"][0]:
                self.assertAlmostEqual(11.0, thickness);
            else:
                self.assertAlmostEqual(0.5, thickness);


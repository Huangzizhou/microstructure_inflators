import unittest
from numpy.linalg import norm
from math import sqrt

from utils.find_file import find_file
from core.WireNetwork import WireNetwork
from parameter.ParameterFactory import ParameterFactory
from ParameterHandler import ParameterHandler

class ParameterHandlerTest(unittest.TestCase):
    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        wire_network = WireNetwork();
        wire_network.load_from_file(wire_file);
        return wire_network;

    def create_parameters(self, wire_network, config):
        factory = ParameterFactory(wire_network);
        factory.create_parameters_from_dict(config);
        return factory.parameters;

    def test_vertex_thickness(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
        config = {
                "thickness": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0],
                    "thickness": [1.0],
                    "default": 0.5
                    }
                }
        parameters = self.create_parameters(wire_network, config);
        handler = ParameterHandler(wire_network);
        handler.convert_to_attributes(parameters);

        self.assertTrue("vertex_thickness" in wire_network.attributes);
        self.assertFalse("edge_thickness" in wire_network.attributes);
        vertex_thickness = wire_network.attributes["vertex_thickness"];
        self.assertEqual(wire_network.num_vertices, len(vertex_thickness));
        self.assertEqual(0.5, min(vertex_thickness));
        self.assertEqual(1.0, max(vertex_thickness));

    def test_edge_thickness(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
        config = {
                "thickness": {
                    "type": "edge_orbit",
                    "effective_orbits": [0],
                    "thickness": [1.0],
                    "default": 0.5
                    }
                }
        parameters = self.create_parameters(wire_network, config);
        handler = ParameterHandler(wire_network);
        handler.convert_to_attributes(parameters);

        self.assertFalse("vertex_thickness" in wire_network.attributes);
        self.assertTrue("edge_thickness" in wire_network.attributes);
        edge_thickness = wire_network.attributes["edge_thickness"];
        self.assertEqual(wire_network.num_edges, len(edge_thickness));
        self.assertEqual(0.5, min(edge_thickness));
        self.assertEqual(1.0, max(edge_thickness));

    def test_vertex_offset(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
        config = {
                "vertex_offset": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0],
                    "offset_percentages": [[1.0, 1.0, 1.0]]
                    }
                }
        parameters = self.create_parameters(wire_network, config);
        handler = ParameterHandler(wire_network);
        handler.convert_to_attributes(parameters);

        self.assertTrue("vertex_offset" in wire_network.attributes);
        offsets = wire_network.attributes["vertex_offset"];
        self.assertEqual(wire_network.num_vertices, len(offsets));
        offset_length = norm(offsets, axis=1);
        self.assertEqual(0.0, min(offset_length));
        self.assertAlmostEqual(sqrt(3.0), max(offset_length));


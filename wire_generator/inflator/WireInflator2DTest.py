import unittest
import numpy as np

from core.WireNetwork import WireNetwork
from parameter.ParameterFactory import ParameterFactory
from utils.find_file import find_file

from WireInflator2D import WireInflator2D

class WireInflator2DTest(unittest.TestCase):
    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        wire_network = WireNetwork();
        wire_network.load_from_file(wire_file);
        return wire_network;

    def create_parameters(self, wire_network, config):
        factory = ParameterFactory(wire_network);
        factory.create_parameters_from_dict(config);
        return factory.parameters;

    def create_config(self):
        config = {
                "thickness": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0],
                    "thickness": [1.0],
                    "default": 0.5
                    },
                "vertex_offset": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0],
                    "offset_percentages": [[0.0, 0.0]]
                    }
                }
        return config;

    def test_periodic(self):
        wire_network = self.load_wire("patterns/2D/box.wire");
        config = self.create_config();
        parameters = self.create_parameters(wire_network, config);
        inflator = WireInflator2D(wire_network, parameters);
        inflator.tile_periodic(np.zeros(2), np.ones(2) * 5);
        mesh = inflator.mesh;
        self.assertGreater(mesh.get_num_vertices(), 0);
        self.assertGreater(mesh.get_num_faces(), 0);

    def test_box(self):
        wire_network = self.load_wire("patterns/2D/box.wire");
        config = self.create_config();
        parameters = self.create_parameters(wire_network, config);
        inflator = WireInflator2D(wire_network, parameters);
        inflator.tile_box(np.zeros(2), np.ones(2) * 10, 2, 2);
        mesh = inflator.mesh;
        self.assertGreater(mesh.get_num_vertices(), 0);
        self.assertGreater(mesh.get_num_faces(), 0);

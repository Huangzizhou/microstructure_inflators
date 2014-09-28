import unittest
import numpy as np
from utils.find_file import find_file
from ParameterFactory import ParameterFactory
from core.WireNetwork import WireNetwork
from VertexThicknessParameter import VertexThicknessParameter
from EdgeThicknessParameter import EdgeThicknessParameter
from VertexOffsetParameter import VertexOffsetParameter

class ParameterFactoryTest(unittest.TestCase):
    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        wire_network = WireNetwork();
        wire_network.load_from_file(wire_file);
        return wire_network;

    def test_creation(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        param_factory = ParameterFactory(wire_network);

    def test_init_with_config(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
        param_factory = ParameterFactory(wire_network, 0.5);
        config = {
                "vertex_offset": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0,1],
                    "offset_percentages": [[0.1, 0.2, 0.3], ["{x}", "{y}", "{x}"]],
                    },
                "thickness": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0, 1],
                    "thickness": [0.1, "{x}"],
                    "default": 0.5
                    }
                }
        param_factory.create_parameters_from_dict(config);
        parameters = param_factory.parameters;

        vertex_orbit = wire_network.attributes["symmetry_vertex_orbit"];
        num_orbits = len(np.unique(vertex_orbit));
        self.assertEqual(num_orbits*2, len(parameters));

        for i,param in enumerate(parameters):
            if i<num_orbits:
                self.assertIsInstance(param, VertexThicknessParameter);
                value = param.evaluate(x=1.0, y=2.0).tolist();
                if i == 0:
                    self.assertEqual([0.1], value);
                elif i == 1:
                    self.assertEqual([1.0], value);
                else:
                    self.assertEqual([0.5], value);
            else:
                self.assertIsInstance(param, VertexOffsetParameter);
                value = param.evaluate(x=1.0, y=2.0).tolist();
                if i%num_orbits == 0:
                    self.assertEqual([0.1, 0.2, 0.3], value);
                elif i%num_orbits == 1:
                    self.assertEqual([1.0, 2.0, 1.0], value);
                else:
                    self.assertEqual([0.0, 0.0, 0.0], value);


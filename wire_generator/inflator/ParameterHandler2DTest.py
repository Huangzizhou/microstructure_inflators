import unittest
import numpy as np
import os

import core.PyWireInflator2DSetting
from core.WireNetwork import WireNetwork
from parameter.ParameterFactory import ParameterFactory
from utils.find_file import find_file
from wire_io.WireWriter import WireWriter

import PyWireInflator2D
from ParameterHandler2D import ParameterHandler2D

class ParameterHandler2DTest(unittest.TestCase):
    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        wire_network = WireNetwork();
        wire_network.load_from_file(wire_file);
        return wire_network;

    def load_inflator(self, wire_network):
        self.assertTrue(2, wire_network.dim);
        tmp_wire_file = "/tmp/test.wire";

        vertices = wire_network.vertices;
        edges = wire_network.edges;
        vertices = np.hstack((vertices, np.zeros((
            wire_network.num_vertices, 1))));

        writer = WireWriter(tmp_wire_file);
        writer.write(vertices, edges);

        inflator = PyWireInflator2D.WireInflatorFacade(tmp_wire_file);
        os.remove(tmp_wire_file);
        return inflator;

    def create_parameters(self, wire_network, config):
        factory = ParameterFactory(wire_network);
        factory.create_parameters_from_dict(config);
        return factory.parameters;

    def test_const_thickness_parameter(self):
        wire_network = self.load_wire("patterns/2D/box.wire");
        inflator = self.load_inflator(wire_network);
        config = {
                "thickness": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0],
                    "thickness": [1.0],
                    "default": 10.5
                    }
                }

        parameters = self.create_parameters(wire_network, config);
        handler = ParameterHandler2D(wire_network, inflator);
        param = handler.convert_to_flattened_parameters(parameters);

        num_parameters = inflator.get_num_parameters();
        self.assertEqual(num_parameters, len(param));

        vertex_orbits = wire_network.attributes["orthotropic_symmetry_vertex_orbit"];
        for i in range(num_parameters):
            affected_vertices = inflator.get_affected_vertex_orbit(i).ravel();
            orbit_ids = vertex_orbits[affected_vertices];
            self.assertTrue(np.all(orbit_ids == orbit_ids[0]));
            orbit_id = orbit_ids[0];
            affected_vertices_2 = np.arange(wire_network.num_vertices,
                    dtype=int)[vertex_orbits == orbit_id];
            self.assertSetEqual(set(affected_vertices), set(affected_vertices_2));

            param_type = inflator.get_parameter_type(i);
            if param_type == PyWireInflator2D.WireInflatorFacade.THICKNESS:
                if orbit_id == 0:
                    self.assertEqual(1.0, param[i]);
                else:
                    self.assertEqual(10.5, param[i]);
            elif param_type == PyWireInflator2D.WireInflatorFacade.VERTEX_OFFSET:
                self.assertEqual(0.0, param[i]);
            else:
                self.fail("Unknown parameter type: {}".format(param_type));

    def test_const_offset(self):
        wire_network = self.load_wire("patterns/2D/box.wire");
        inflator = self.load_inflator(wire_network);
        config = {
                "thickness": {
                    "type": "vertex_orbit",
                    "effective_orbits": [],
                    "thickness": [],
                    "default": 0.5
                    },
                "vertex_offset": {
                    "type": "vertex_orbit",
                    "effective_orbits": [0],
                    "offset_percentages": [[1.0, 1.0]]
                    }
                }

        parameters = self.create_parameters(wire_network, config);
        handler = ParameterHandler2D(wire_network, inflator);
        param = handler.convert_to_flattened_parameters(parameters);

        num_parameters = inflator.get_num_parameters();
        self.assertEqual(num_parameters, len(param));

        vertex_orbits = wire_network.attributes["orthotropic_symmetry_vertex_orbit"];
        for i in range(num_parameters):
            affected_vertices = inflator.get_affected_vertex_orbit(i).ravel();
            orbit_ids = vertex_orbits[affected_vertices];
            self.assertTrue(np.all(orbit_ids == orbit_ids[0]));
            orbit_id = orbit_ids[0];
            param_type = inflator.get_parameter_type(i);
            if param_type == PyWireInflator2D.WireInflatorFacade.VERTEX_OFFSET:
                if orbit_id == 0:
                    self.assertEqual(1.0, param[i]);
                else:
                    self.assertEqual(0.0, param[i]);
            elif param_type == PyWireInflator2D.WireInflatorFacade.THICKNESS:
                self.assertEqual(0.5, param[i]);
            else:
                self.fail("Unknown parameter type: {}".format(param_type));


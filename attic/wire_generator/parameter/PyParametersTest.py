import unittest
import json
import numpy as np
from numpy.linalg import norm
from utils.find_file import find_file
from core.WireNetwork import WireNetwork
import os
from PyParameters import PyParameters

class PyParametersTest(unittest.TestCase):
    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        wire_network = WireNetwork();
        wire_network.load_from_file(wire_file);
        return wire_network;

    def test_empty(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameters = PyParameters(wire_network, 0.5);
        self.assertEqual(0, parameters.num_dofs);
        self.assertEqual(0, parameters.num_thickness_dofs);
        self.assertEqual(0, parameters.num_offset_dofs);

    def test_modifier_file(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
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
        parameters = PyParameters(wire_network, 0.5);
        parameters.load_modifier(config);

        self.assertEqual(6, parameters.num_dofs);
        self.assertEqual(2, parameters.num_thickness_dofs);
        self.assertEqual(4, parameters.num_offset_dofs);

        self.assertAlmostEqual(0.1, parameters.dofs[0]);
        self.assertAlmostEqual(0.0, parameters.dofs[1]);
        self.assertAlmostEqual(0.1, parameters.dofs[2]);
        self.assertAlmostEqual(0.2, parameters.dofs[3]);
        self.assertAlmostEqual(0.3, parameters.dofs[4]);
        self.assertAlmostEqual(0.0, parameters.dofs[5]);

    def test_dof_file(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
        dof_setting = {
                "dof_type": "isotropic",
                "thickness_type": "vertex",
                "dof": [0.3, 0.4, 0.5, 0.1, 0.2]
                };
        dof_file = "/tmp/tmp.dof";
        with open(dof_file, 'w') as fout:
            json.dump(dof_setting, fout);
        parameters = PyParameters(wire_network, 0.5);
        parameters.load_dof_file(dof_file);

        self.assertEqual(5, parameters.num_dofs);
        self.assertEqual(3, parameters.num_thickness_dofs);
        self.assertEqual(2, parameters.num_offset_dofs);

        dof_diff = np.array(dof_setting["dof"]) - parameters.dofs;
        self.assertAlmostEqual(0.0, norm(dof_diff));

        os.remove(dof_file);

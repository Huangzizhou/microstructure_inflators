import unittest
import numpy as np

from utils.find_file import find_file
from core.WireNetwork import WireNetwork
from WireSymmetryAttribute import WireSymmetryAttribute

class WireSymmetryAttributeTest(unittest.TestCase):
    def setUp(self):
        self.attribute = WireSymmetryAttribute();

    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        wire_network = WireNetwork();
        wire_network.load_from_file(wire_file);
        return wire_network;

    def assert_symmetry_is_invertible(self, wire_network):
        for symm in self.attribute.symmetries:
            vertices = np.array(map(symm, map(symm, wire_network.vertices)));
            self.assertTrue(np.all(vertices == wire_network.vertices));

    def test_2D(self):
        wire_network = self.load_wire("patterns/2D/box.wire");
        self.attribute.compute(wire_network);
        self.assert_symmetry_is_invertible(wire_network);

    def test_3D(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
        self.attribute.compute(wire_network);
        self.assert_symmetry_is_invertible(wire_network);

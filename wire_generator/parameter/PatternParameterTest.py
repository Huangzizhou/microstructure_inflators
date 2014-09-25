import unittest
from utils.find_file import find_file
from core.WireNetwork import WireNetwork
from PatternParameter import PatternParameter

class PatternParameterTest(unittest.TestCase):
    def load_wire(self, wire_file):
        wire_file = find_file(wire_file);
        wire_network = WireNetwork();
        wire_network.load_from_file(wire_file);
        return wire_network;

    def test_formula(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        formula = "x + y"
        parameter = PatternParameter(wire_network);
        parameter.set_formula(formula);

        self.assertEqual(formula, parameter.formula);

import unittest
from WireNetwork import WireNetwork
from WireSymmetryEdgeOrbitAttribute import WireSymmetryEdgeOrbitAttribute

class WireSymmetryEdgeOrbitAttributeTest(unittest.TestCase):
    def load_wire(self, wire_file):
        self.wire_network = WireNetwork();
        self.wire_network.load_from_file(wire_file);
        self.wire_network.attributes.add("symmetry_vertex_orbit");
        self.wire_network.attributes.add("symmetry_edge_orbit");

    def get_orbits(self):
        return self.wire_network.attributes["symmetry_edge_orbit"];

    def test_creation(self):
        self.load_wire("examples/cube.wire");
        self.assertTrue("symmetry_edge_orbit" in self.wire_network.attributes);

    def test_cube(self):
        self.load_wire("examples/cube.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.edges));
        self.assertEqual(3, num_orbits);

    def test_brick5(self):
        self.load_wire("examples/example2.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.edges));
        self.assertEqual(6, num_orbits);

    def test_tet(self):
        self.load_wire("examples/tet.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.edges));
        self.assertEqual(6, num_orbits);

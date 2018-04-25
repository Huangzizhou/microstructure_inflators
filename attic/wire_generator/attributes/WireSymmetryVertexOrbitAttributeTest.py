import unittest
from core.WireNetwork import WireNetwork
from WireSymmetryVertexOrbitAttribute import WireSymmetryVertexOrbitAttribute

class WireSymmetryVertexOrbitAttributeTest(unittest.TestCase):
    def load_wire(self, wire_file):
        self.wire_network = WireNetwork();
        self.wire_network.load_from_file(wire_file);
        self.wire_network.attributes.add("orthotropic_symmetry_vertex_orbit");

    def get_orbits(self):
        return self.wire_network.attributes["orthotropic_symmetry_vertex_orbit"];

    def test_creation(self):
        self.load_wire("examples/cube.wire");
        self.assertTrue("orthotropic_symmetry_vertex_orbit" in self.wire_network.attributes);

    def test_cube(self):
        self.load_wire("examples/cube.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.vertices));
        self.assertEqual(1, num_orbits);

    def test_brick5(self):
        self.load_wire("examples/example2.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.vertices));
        self.assertEqual(7, num_orbits);

    def test_tet(self):
        self.load_wire("examples/tet.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.vertices));
        self.assertEqual(1, num_orbits);

    def test_cube_isotropic(self):
        self.load_wire("examples/cube.wire");
        self.wire_network.attributes.add("isotropic_symmetry_vertex_orbit");
        orbits = self.wire_network.attributes[
                "isotropic_symmetry_vertex_orbit"];

        num_orbits = len(set(orbits));
        self.assertEqual(len(orbits), len(self.wire_network.vertices));
        self.assertEqual(1, num_orbits);

    def test_brick5_isotropic(self):
        self.load_wire("examples/example2.wire");
        self.wire_network.attributes.add("isotropic_symmetry_vertex_orbit");
        orbits = self.wire_network.attributes[
                "isotropic_symmetry_vertex_orbit"];

        num_orbits = len(set(orbits));
        self.assertEqual(len(orbits), len(self.wire_network.vertices));
        self.assertEqual(3, num_orbits);
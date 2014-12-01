import unittest
from core.WireNetwork import WireNetwork
from WireSymmetryEdgeOrbitAttribute import WireSymmetryEdgeOrbitAttribute

import numpy as np
from numpy.linalg import norm

class WireSymmetryEdgeOrbitAttributeTest(unittest.TestCase):
    def load_wire(self, wire_file):
        self.wire_network = WireNetwork();
        self.wire_network.load_from_file(wire_file);
        self.wire_network.attributes.add("orthotropic_symmetry_vertex_orbit");
        self.wire_network.attributes.add("orthotropic_symmetry_edge_orbit");

    def get_orbits(self):
        return self.wire_network.attributes["orthotropic_symmetry_edge_orbit"];

    def assertSymmetricVectors(self, vectors):
        for i in range(len(vectors)-1):
            vi = vectors[i];
            vj = vectors[i+1];
            diff = np.absolute(vi) - np.absolute(vj);
            self.assertAlmostEqual(0.0, norm(diff));

    def assertOrbitIsValid(self, orbit):
        edges = self.wire_network.edges[orbit];
        edge_vectors = self.wire_network.vertices[edges[:,0]] -\
                self.wire_network.vertices[edges[:,1]];
        self.assertSymmetricVectors(edge_vectors);

    def assertOrbitsAreValid(self, orbits):
        orbit_map = {};
        for i,oi in enumerate(orbits):
            orbit_map[oi] = orbit_map.get(oi, []) + [i];

        for i,orbit in orbit_map.iteritems():
            self.assertOrbitIsValid(orbit);

    def test_creation(self):
        self.load_wire("examples/cube.wire");
        self.assertTrue("orthotropic_symmetry_edge_orbit" in self.wire_network.attributes);

    def test_cube(self):
        self.load_wire("examples/cube.wire");
        orbits = self.get_orbits();

        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.edges));
        self.assertEqual(3, num_orbits);
        self.assertOrbitsAreValid(orbits);


    def test_brick5(self):
        self.load_wire("examples/example2.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.edges));
        self.assertEqual(6, num_orbits);
        self.assertOrbitsAreValid(orbits);

    def test_tet(self):
        self.load_wire("examples/tet.wire");
        orbits = self.get_orbits();
        num_orbits = len(set(orbits));

        self.assertEqual(len(orbits), len(self.wire_network.edges));
        self.assertEqual(6, num_orbits);
        self.assertOrbitsAreValid(orbits);

    def test_brick5_isotropic(self):
        self.load_wire("examples/example2.wire");
        self.wire_network.attributes.add("isotropic_symmetry_edge_orbit");
        orbits = self.wire_network.attributes["isotropic_symmetry_edge_orbit"];

        num_orbits = len(set(orbits));
        self.assertEqual(len(orbits), len(self.wire_network.edges));
        self.assertEqual(2, num_orbits);

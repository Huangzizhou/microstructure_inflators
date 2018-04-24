from ThicknessParameterTest import ThicknessParameterTest
from EdgeThicknessParameter import EdgeThicknessParameter

class EdgeThicknessParameterTest(ThicknessParameterTest):
    def names_test(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameter = EdgeThicknessParameter(wire_network, 0, 1.0);
        names = parameter.names;
        self.assertEqual(["edge_orbit_0_thickness"], names);

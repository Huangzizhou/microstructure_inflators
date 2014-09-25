from ThicknessParameterTest import ThicknessParameterTest
from VertexThicknessParameter import VertexThicknessParameter

class VertexThicknessParameterTest(ThicknessParameterTest):
    def names_test(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameter = VertexThicknessParameter(wire_network, 1, 1.0);
        names = parameter.names;
        self.assertEqual(["vertex_orbit_1_thickness"], names);

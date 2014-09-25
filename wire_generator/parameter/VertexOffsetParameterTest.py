import unittest
import numpy as np
from PatternParameterTest import PatternParameterTest
from VertexOffsetParameter import VertexOffsetParameter

class VertexOffsetParameterTest(PatternParameterTest):
    def test_default_3D(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        self.assertEqual(3, len(parameter.default_offset));

    @unittest.skip("debug")
    def test_default_2D(self):
        wire_network = self.load_wire("patterns/2D/box.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        self.assertEqual(2, len(parameter.default_offset));

    def test_evaluate(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        formula = ["{x}+{y}", "{x}-{y}"];
        parameter.set_formula(formula);
        result = parameter.evaluate(x=1, y=2);
        self.assertEqual([3.0, -1.0], result.tolist());

    def test_evaluate_default(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        result = parameter.evaluate();
        self.assertEqual([0.0, 0.0, 0.0], result.tolist());

    def test_evaulate_fixed_offset(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        formula = [3.0, -1.0];
        parameter.set_formula(formula);
        result = parameter.evaluate();
        self.assertEqual([3.0, -1.0], result.tolist());

    def test_names(self):
        wire_network = self.load_wire("patterns/3D/star.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        names = parameter.names;
        self.assertEqual([
            "vertex_orbit_0_offset_x",
            "vertex_orbit_0_offset_y",
            "vertex_orbit_0_offset_z" ], names);

    def test_dof_mask(self):
        wire_network = self.load_wire("patterns/3D/cube.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        mask = parameter.dof_mask;
        self.assertTrue(np.all(mask));

    def test_dof_mask(self):
        wire_network = self.load_wire("patterns/3D/brick5.wire");
        parameter = VertexOffsetParameter(wire_network, 0);
        mask = parameter.dof_mask;
        self.assertFalse(np.all(mask));


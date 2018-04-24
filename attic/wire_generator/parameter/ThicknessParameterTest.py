from PatternParameterTest import PatternParameterTest
from ThicknessParameter import ThicknessParameter

class ThicknessParameterTest(PatternParameterTest):
    def test_evaluate(self):
        wire_network = self.load_wire("patterns/3D/star.wire");

        formula = "{x}+{y}"
        parameter = ThicknessParameter(wire_network, 0, 1.0);
        parameter.set_formula(formula);
        self.assertEqual(0, parameter.orbit_id);
        self.assertEqual(1.0, parameter.default_thickness);
        self.assertEqual(formula, parameter.formula);

        result = parameter.evaluate(x=1, y=2);
        self.assertEqual(1, len(result));
        self.assertEqual(3, result[0]);

    def test_evaluate_default(self):
        wire_network = self.load_wire("patterns/3D/star.wire");

        parameter = ThicknessParameter(wire_network, 0, 1.0);
        result = parameter.evaluate();
        self.assertEqual([1.0], result.tolist());

    def test_evaluate_fixed_thickness(self):
        wire_network = self.load_wire("patterns/3D/star.wire");

        parameter = ThicknessParameter(wire_network, 0, 1.0);
        parameter.set_formula(0.5);
        result = parameter.evaluate();
        self.assertEqual([0.5], result.tolist());

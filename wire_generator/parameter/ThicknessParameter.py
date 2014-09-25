import numpy as np
from PatternParameter import PatternParameter

class ThicknessParameter(PatternParameter):
    def __init__(self, wire_network, orbit_id, default_thickness = 0.5):
        super(ThicknessParameter, self).__init__(wire_network);
        self.default_thickness = default_thickness;
        self.orbit_id = orbit_id;

    def evaluate(self, **kwargs):
        if hasattr(self, "formula"):
            if isinstance(self.formula, (str, unicode)):
                formula = self.formula.format(**kwargs);
                return np.array([eval(formula)]);
            elif isinstance(self.formula, float):
                return np.array([self.formula]);
            else:
                raise RuntimeError("Unable to evaluate formula {}".format(
                    self.formula));
        else:
            return np.array([self.default_thickness]);



import numpy as np

class PatternParameter(object):
    def __init__(self, wire_network):
        self.wire_network = wire_network;

    def evaluate(self, **kwargs):
        raise NotImplementedError("This method is abstract.");

    def set_formula(self, formula):
        self.formula = formula;

    @property
    def names(self):
        raise NotImplementedError("This method is abstract.");

    @property
    def dof_mask(self):
        return np.zeros(len(self.names), dtype=True);

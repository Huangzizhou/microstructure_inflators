class OptimizationSetting(object):
    def evaluate_objective(self, parameters):
        raise NotImplementedError("This method is abstract");

    def evaluate_gradient(self, parameters):
        raise NotImplementedError("This method is abstract");

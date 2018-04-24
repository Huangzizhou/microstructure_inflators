class PolynomialFitter(object):
    @classmethod
    def create(cls, poly_type):
        if (poly_type == "quadratic"):
            from QuadraticFitter import QuadraticFitter
            return QuadraticFitter();
        elif (poly_type == "separable_quadratic"):
            from SeparableQuadraticFitter import SeparableQuadraticFitter
            return SeparableQuadraticFitter();
        else:
            raise NotImplementedError("Unsupported fitter type: {}".format(poly_type));

    def fit(self, x, y):
        raise NotImplementedError("This method is abstract");

    def evaluate(self, x):
        raise NotImplementedError("This method is abstract");

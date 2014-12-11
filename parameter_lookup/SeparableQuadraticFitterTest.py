import unittest
from PolynomialFitter import PolynomialFitter

import numpy as np
from numpy.linalg import norm

class SeparableQuadraticFitterTest(unittest.TestCase):
    def setUp(self):
        self.x = np.random.rand(10, 2);
        self.y = np.arange(10) * 10;
        #self.y = np.multiply(np.square(self.x[:, 0]), self.x[:, 1]);
        self.y = self.y.reshape((-1, 1));
        for x,y in zip(self.x, self.y):
            print(x, y);

    def test_separable_quadratic_fit(self):
        fitter = PolynomialFitter.create("separable_quadratic");
        fitter.fit(self.x, self.y);
        fitted_y = fitter.evaluate(self.x);

        err = norm(fitted_y - self.y);
        self.assertAlmostEqual(0.0, err);

    def test_evaluate(self):
        fitter = PolynomialFitter.create("separable_quadratic");
        fitter.x_num_dof = 2;
        fitter.y_num_dof = 1;
        fitter.sol = np.ones((3, self.x.shape[1], self.y.shape[1]));
        y = fitter.evaluate(self.x);

import unittest
import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize

class ScipyOptimizeTest(unittest.TestCase):
    def optimize(self, f, init_x, method="L-BFGS-B", bounds=None, disp=False):
        result = minimize(f, init_x,
                method=method, jac=True,
                bounds=bounds,
                options={
                    "maxiter": 10,
                    "disp": disp
                    });
        return result;

    def test_easy(self):
        def f(x):
            return norm(x)**2, 2*x;

        result = self.optimize(f, [0.5, 0.5], bounds=[(-1.0, 1.0), (-1.0, 1.0)]);
        self.assertEqual(0.0, norm(result.x));

    def test_easy_2(self):
        def f(x):
            return norm(x), 2*x;

        result = self.optimize(f, [0.5, 0.5], bounds=[(0.1, 1.0), (0.1, 1.0)]);
        self.assertEqual(0.1, np.max(result.x));
        self.assertEqual(0.1, np.min(result.x));

    def test_easy_3(self):
        def f(x):
            return x[0], np.array([1.0, 0.0]);

        result = self.optimize(f, [0.5, 0.5], bounds=[(0.1, 1.0), (0.1, 1.0)]);
        self.assertEqual(0.1, result.x[0]);

    def test_easy_4(self):
        def f(x):
            return x[0]**5, np.array([5*x[0]**4, 0.0]);

        result = self.optimize(f, [0.5, 0.5], bounds=[(0.1, 1.0), (0.1, 1.0)],
                disp=False);
        self.assertEqual(0.1, result.x[0]);

    def test_concave_down(self):
        def f(x):
            return -norm(x)**2, -2*x;

        result = self.optimize(f, [0.5, 0.5], bounds=[(0.1, 1.0), (0.1, 1.0)]);
        self.assertEqual(1.0, np.max(result.x));
        self.assertEqual(1.0, np.min(result.x));


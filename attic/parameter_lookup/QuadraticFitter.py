import numpy as np
from numpy.linalg import lstsq
from PolynomialFitter import PolynomialFitter

class QuadraticFitter(PolynomialFitter):
    def fit(self, x, y):
        n = x.shape[0];
        assert(n == y.shape[0]);

        self.x_dof = x.shape[1];
        self.y_dof = y.shape[1];

        coeff_mat = [np.ones((n, 1))];
        coeff_mat.append(self.__linear_coeffs(x));
        coeff_mat.append(self.__quadratic_coeffs(x));
        #coeff_mat.append(self.__quartic_coeffs(x));

        coeff_mat = np.hstack(coeff_mat);
        sol, residual, rank, singular_vals = lstsq(coeff_mat, y);

        self.sol= sol;

    def evaluate(self, x):
        x = x.reshape((-1, self.x_dof), order="C");
        n = x.shape[0];

        coeff_mat = [np.ones((n, 1))];
        coeff_mat.append(self.__linear_coeffs(x));
        coeff_mat.append(self.__quadratic_coeffs(x));
        #coeff_mat.append(self.__quartic_coeffs(x));

        coeff_mat = np.hstack(coeff_mat);

        return np.dot(coeff_mat, self.sol);

    def __linear_coeffs(self, x):
        return x;

    def __quadratic_coeffs(self, x):
        x_dof = x.shape[1];
        coeffs = [];

        for row in x:
            coeffs.append(np.outer(row, row).reshape((1,-1)));

        return np.vstack(coeffs);

    def __quartic_coeffs(self, x):
        x_dof = x.shape[1];
        coeffs = [];
        for row in x:
            quadratic_terms = np.outer(row, row).ravel();
            quartic_terms = np.outer(quadratic_terms, quadratic_terms).ravel();
            coeffs.append(quartic_terms);
        return np.vstack(coeffs);

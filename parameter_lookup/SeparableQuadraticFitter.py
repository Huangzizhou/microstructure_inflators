import numpy as np
from numpy.linalg import lstsq, norm
from PolynomialFitter import PolynomialFitter

class SeparableQuadraticFitter(PolynomialFitter):
    def fit_old(self, x, y):
        n = x.shape[0];
        self.x_num_dof = x.shape[1];
        self.y_num_dof = y.shape[1];

        self.sol = np.ones((3, self.x_num_dof, self.y_num_dof));
        max_itr = 100;
        for i in range(max_itr):
            sol = np.copy(self.sol);
            for j in range(self.x_num_dof):
                self.__fit_some_variables(j, x, y);

            change = np.amax(np.absolute(sol - self.sol));
            err = np.amax(np.absolute(self.evaluate(x) - y));

            if change < 1e-6 or err < 1e-6:
                break;

    def fit(self, x, y):
        n = x.shape[0];
        self.x_num_dof = x.shape[1];
        self.y_num_dof = y.shape[1];

        self.sol = np.ones((3, self.x_num_dof, self.y_num_dof));
        from scipy.optimize import leastsq
        sol = leastsq(self.compute_residual, self.sol.ravel(order="C"),
                args=(x, y))[0];
        self.sol = sol.reshape((3, self.x_num_dof, self.y_num_dof), order="C");

    def evaluate(self, x, x_mask=None):
        x = x.reshape((-1, self.x_num_dof));
        const_terms = np.ones_like(x);
        linear_terms = x;
        quadratic_terms = np.square(x);

        if x_mask is None:
            x_mask = np.ones(self.x_num_dof, dtype=bool);

        results = [];
        for i in range(self.y_num_dof):
            sol = self.sol[:, :, i];
            r  = const_terms * sol[0, :].ravel();
            r += linear_terms * sol[1, :].ravel();
            r += quadratic_terms * sol[2, :].ravel();
            r = r[:, x_mask];
            r = np.prod(r, axis=1).reshape((-1, 1));
            results.append(r);
        results = np.hstack(results);
        return results;

    def compute_residual(self, sol, x, y):
        self.sol = sol.reshape((3, self.x_num_dof, self.y_num_dof), order="C");
        res = (self.evaluate(x) - y).ravel(order="C");
        return res;

    def __fit_some_variables(self, dof_idx, x, y):
        num_dof = x.shape[1];
        unknown = x[:, dof_idx].reshape((-1, 1));

        unknown_const_terms = np.ones_like(unknown);
        unknown_linear_terms = unknown;
        unknown_quadratic_terms = np.square(unknown);
        unknown_coeff_mat = np.hstack([
            unknown_const_terms,
            unknown_linear_terms,
            unknown_quadratic_terms ]);

        var_mask = np.arange(num_dof, dtype=int) != dof_idx;
        known_lhs = self.evaluate(x, var_mask);
        #rhs = np.divide(y, known_lhs);

        err_before = self.__compute_error(x, y);

        for i in range(self.y_num_dof):
            A = np.copy(unknown_coeff_mat);
            A *= known_lhs[:,i].reshape((-1, 1));
            sol_before = np.copy(self.sol[:, dof_idx, i]);
            residual_before = norm(np.dot(A, sol_before) - y[:,i])**2;
            sol, residual, rank, singular_vals = lstsq(A, y[:,i]);
            self.sol[:, dof_idx, i] = sol.ravel();
            residual_after = norm(np.dot(A, self.sol[:, dof_idx, i]) - y[:,i])**2;
            if (residual_before < residual_after):
                print(residual_before, residual_after)
                print("========");
                print(sol_before);
                print("========");
                print(sol.ravel());
                print("========");
                assert(False);

            err_after = self.__compute_error(x, y);
            print(err_before, err_after);
            assert(err_before >= err_after);
            err_before = err_after;

        #sol, residual, rank, singular_vals = lstsq(unknown_coeff_mat, rhs);
        #self.sol[:, dof_idx, :] = sol;

        #err_after = self.__compute_error(x, y);
        #print(err_before, err_after);
        #assert(err_before >= err_after);

    def __compute_error(self, x, y):
        fitted_y = self.evaluate(x);
        return norm(fitted_y - y);

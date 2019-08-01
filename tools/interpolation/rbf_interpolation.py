import math

import numpy as np
import scipy
from scipy import optimize

class RBFInterpolation:

    def __init__(self, x1 = [], x2 = [], p = [], d1 = 5, d2 = 5, epsilon = 0.5, coeffs = []):
        self.d1 = d1
        self.d2 = d2

        self.min1 = -1.0
        self.max1 =  1.0
        self.min2 = -1.0
        self.max2 = 1.0

        self.epsilon = epsilon

        self.dt1 = (self.max1 - self.min1) / (self.d1 - 1)
        self.dt2 = (self.max2 - self.min2) / (self.d2 - 1)

        if len(coeffs) > 0:
            self.coeffs = coeffs
            return

        # create rhs vector
        B = np.array(p)

        # populating A
        A = np.zeros((len(x1), self.d1 * self.d2))

        x1 = np.array(x1)
        x2 = np.array(x2)
        for i in range(0, self.d1):
            for j in range(0, self.d2):
                # Find rbf point correspond to i and j
                xi = self.min1 + i * self.dt1
                xj = self.min2 + j * self.dt1

                r = np.sqrt(np.power(x1 - xi, 2) + np.power(x2 - xj, 2))

                A[:, self.d2 * i + j] = np.exp(- np.power((self.epsilon * r), 2))


        # solving least squares
        x = scipy.optimize.lsq_linear(A, B)

        # saving coefficients
        self.coeffs = x.x


    def __call__(self, x1, x2):
        results = np.zeros(len(x1))

        coeff_idx = 0
        for i in range(0, self.d1):
            for j in range(0, self.d2):
                # Find rbf point correspond to i and j
                xi = self.min1 + i * self.dt1
                xj = self.min2 + j * self.dt1

                r = np.sqrt(np.power(x1 - xi, 2) + np.power(x2 - xj, 2))

                results += self.coeffs[coeff_idx] * np.exp(- np.power((self.epsilon * r), 2))

                coeff_idx += 1
        return results
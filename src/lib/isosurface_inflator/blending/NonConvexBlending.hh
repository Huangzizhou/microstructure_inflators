//
// Created by Davi Colli Tozoni on 8/11/18.
//

#ifndef NONCONVEXBLENDING_H
#define NONCONVEXBLENDING_H

namespace SignedDistance {
    // This function corresponds to an approximation of the absolute function. The nonconvex name is refers to the polynomial
    // that gives rise to this final equation (check documentation in doc/blending.pdf). Basically, this approximation does
    // not guarantee same sign curvature on the whole blending region.
    // Function was automatically extracted from Maple script doc/KS_nonconvex_8parameters.mw
    template<typename Real>
    Real
    compute_nonconvex_fake_log_2cosh(Real &x, Real &p, const Real &a0, const Real &a2, const Real &a4, const Real &a6, const Real &a8, const Real &a10, const Real &a12,
                                     Real a14) {
        return ((-0.40040e5 * a14 * pow(p, 0.18e2) * pow(x, 0.18e2) +
                 0.204204e6 * (a14 - a12 / 0.4e1) * pow(x, 0.16e2) * pow(p, 0.16e2) -
                 0.67320e5 * pow(x, 0.14e2) * (a10 - 0.4e1 * a12) * pow(p, 0.14e2) -
                 0.92820e5 * pow(x, 0.12e2) * (a8 - 0.4e1 * a10) * pow(p, 0.12e2) -
                 0.136136e6 * pow(x, 0.10e2) * (a6 - 0.4e1 * a8) * pow(p, 0.10e2) -
                 0.218790e6 * pow(x, 0.8e1) * (a4 - 0.4e1 * a6) * pow(p, 0.8e1) -
                 0.408408e6 * pow(x, 0.6e1) * (a2 - 0.4e1 * a4) * pow(p, 0.6e1) -
                 0.1021020e7 * pow(x, 0.4e1) * (a0 - 0.4e1 * a2) * pow(p, 0.4e1) + 0.24504480e8 * a0 * p * p * x * x +
                 0.22304522240e11 * a14 + 0.49008960e8 * a0 + 0.65345280e8 * a2 + 0.130690560e9 * a4 +
                 0.313657344e9 * a6 + 0.836419584e9 * a8 + 0.2389770240e10 * a10 + 0.7169310720e10 * a12) /
                (0.49201152e8 * a14 + 0.16084992e8 * a12 + 0.5483520e7 * a10 + 0.1980160e7 * a8 + 0.777920e6 * a6 +
                 0.350064e6 * a4 + 0.204204e6 * a2 + 0.255255e6 * a0) / p / 0.256e3);
    }

    // same as function presented above, but receives vector of coefficients.
    template<typename Real>
    Real compute_nonconvex_fake_log_2cosh(Real x, Real p, const std::vector<Real> &poly_coeffs) {
        Real result = 0.0;

        // Building default values. To create the data structure needed for automatic differentiation, the easiest way is
        // to multiply the default values of coefficients by 1, already containing the derivatives data structure.
        Real one = p / p;
        Real da0 = 0.01 * one;
        Real da2 = 0.0001 * one;
        Real da4 = 0.0 * one;
        Real da6 = 0.0 * one;
        Real da8 = 0.0 * one;
        Real da10 = 0.0 * one;
        Real da12 = 0.0 * one;
        Real da14 = 0.0 * one;

        // If x is outside the blending range, we assume the value of the absolute function
        if (x < -2.0 / p)
            return -x;
        if (x > 2.0 / p)
            return x;

        // Because, aparently, the version I have of Eigen3 does not take care of operations where one of the variables is
        // not constant and not of the automatic diff type, and, also, because I extracted the functions directly from maple,
        // I thought it would be clearer (mainly in numbering) to keep parameters as individual variables instead of an array.
        switch (poly_coeffs.size()) {
            case 0:
                result = compute_nonconvex_fake_log_2cosh(x, p, da0, da2, da4, da6, da8, da10, da12, da14);
                break;
            case 1:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], da2, da4, da6, da8, da10, da12, da14);
                break;
            case 2:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], da4, da6, da8, da10,
                                                          da12, da14);
                break;
            case 3:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2], da6,
                                                          da8, da10, da12, da14);
                break;
            case 4:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2],
                                                          poly_coeffs[3], da8, da10, da12, da14);
                break;
            case 5:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2],
                                                          poly_coeffs[3], poly_coeffs[4], da10, da12, da14);
                break;
            case 6:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2],
                                                          poly_coeffs[3], poly_coeffs[4], poly_coeffs[5], da12, da14);
                break;
            case 7:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2],
                                                          poly_coeffs[3], poly_coeffs[4], poly_coeffs[5],
                                                          poly_coeffs[6], da14);
                break;
            case 8:
                result = compute_nonconvex_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2],
                                                          poly_coeffs[3], poly_coeffs[4], poly_coeffs[5],
                                                          poly_coeffs[6], poly_coeffs[7]);
                break;
            default:
                throw std::runtime_error("Non recognizable number of extra parameters for blending: " +
                                         std::to_string(poly_coeffs.size()));
        }

        if (hasInvalidDerivatives(result) || std::isnan(stripAutoDiff(result)) || std::isinf(stripAutoDiff(result))) {
            std::cerr << "  result: ";
            std::cerr << result << std::endl;

            std::cerr << "Report derivatives: " << std::endl;
            std::cerr << "  result: ";
            reportDerivatives(std::cerr, result);
            std::cerr << std::endl;
            std::cerr << std::endl;
        }

        return result;
    }


    // Smooth version of min function based on nonconvex approximation of absolute function shown above.
    // This function deals with only two entries
    template<typename Real>
    Real exp_smin_nonconvex_params(Real a, Real b, Real s, const std::vector<Real> &poly_coeffs) {
        if (s < 1e-10) {
            return std::min(a, b);
        }

        Real epsilon = (a - b) / 2.0;
        Real p = 1.0 / s;
        Real result = (a + b) / 2.0 - compute_nonconvex_fake_log_2cosh(epsilon, p, poly_coeffs);

        if (hasInvalidDerivatives(result) || std::isnan(stripAutoDiff(result)) || std::isinf(stripAutoDiff(result))) {
            std::cerr << "  result: ";
            std::cerr << result << std::endl;
            std::cerr << "  s: ";
            std::cerr << s << std::endl;
            std::cerr << "  p: ";
            std::cerr << p << std::endl;
            std::cerr << "  epsilon: ";
            std::cerr << epsilon << std::endl;

            for (size_t i = 0; i < poly_coeffs.size(); i++) {
                std::cerr << "  poly_coeffs[" << i << "]: ";
                std::cerr << poly_coeffs[i] << std::endl;
            }

            std::cerr << "Report derivatives: " << std::endl;
            std::cerr << "  result: ";
            reportDerivatives(std::cerr, result);
            std::cerr << std::endl;
            std::cerr << std::endl;
            std::cerr << "  s: ";
            reportDerivatives(std::cerr, s);
            std::cerr << std::endl;
            std::cerr << std::endl;
            std::cerr << "  p: ";
            reportDerivatives(std::cerr, p);
            std::cerr << std::endl;
            std::cerr << std::endl;
            std::cerr << "  epsilon: ";
            reportDerivatives(std::cerr, epsilon);
            std::cerr << std::endl;
            std::cerr << std::endl;
        }

        return result;
    }

    // Recursive smooth version of min function based on nonconvex approximation of absolute function shown above.
    // This function deals with an array of values
    template<typename Real>
    Real exp_smin_nonconvex_params(const std::vector<Real> &values, Real s, const std::vector<Real> &poly_coeffs) {
        Real result = 0.0;
        size_t n = values.size();

        if (s <= 1e-10) {
            Real minVal = values[0];
            for (const Real &val : values) {
                minVal = std::min(minVal, val);
            }

            return minVal;
        }

        // Base case
        if (n == 2) {
            result = exp_smin_nonconvex_params<Real>(values[0], values[1], s, poly_coeffs);
            return result;
        }

        // Copy values to array to be used in recursive call.
        // At each iteration i, we want to compute recursively, removing the ith from the values array.
        std::vector<Real> recursionValues(values);
        Real lastElement = values[n - 1];
        recursionValues.resize(n - 1);
        for (size_t i = 0; i < n; i++) {
            Real currentElement = values[i];

            if (i < (n - 1)) {
                // exchange current position by last value
                recursionValues[i] = lastElement;
            }

            Real recursiveS = s / 10.0;
            Real recursiveMin = exp_smin_nonconvex_params<Real>(recursionValues, recursiveS, poly_coeffs);
            result += exp_smin_nonconvex_params<Real>(currentElement, recursiveMin, s, poly_coeffs);

            if (i < (n - 1)) {
                // put back original element
                recursionValues[i] = currentElement;
            }

            // Check if everyhting is alright!
            if (hasInvalidDerivatives(result) || std::isnan(stripAutoDiff(result)) ||
                std::isinf(stripAutoDiff(result))) {
                std::cerr << "  result: ";
                std::cerr << result << std::endl;
                std::cerr << "  s: ";
                std::cerr << s << std::endl;
                std::cerr << "  recursiveMin: ";
                std::cerr << recursiveMin << std::endl;

                std::cerr << "Report derivatives: " << std::endl;
                std::cerr << "  result: ";
                reportDerivatives(std::cerr, result);
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "  s: ";
                reportDerivatives(std::cerr, s);
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "  recursiveMin: ";
                reportDerivatives(std::cerr, recursiveMin);
                std::cerr << std::endl;
                std::cerr << std::endl;

                throw std::runtime_error("Error while computing smooth min with a polynomial convex function");
            }
        }
        result /= (Real) n;

        return result;
    }

}
#endif //NONCONVEXBLENDING_H

//
// Created by Davi Colli Tozoni on 8/11/18.
//

#ifndef POLYNOMIALCONVEXBLENDING_H
#define POLYNOMIALCONVEXBLENDING_H

namespace SignedDistance {

    // This function corresponds to an approximation of the absolute function. The  name is refers to the fact that this is
    // a polynomial and symmetric function. The way this function was produced is explained in the documentation in doc/blending.pdf).
    // Function was automatically extracted from Maple script doc/KS_monotonic_derivative_8parameters.mw
    template<typename Real>
    Real compute_symmetric_fake_log_2cosh(Real x, Real p, const Real &a0, const Real &a2, const Real &a4, const Real &a6, const Real &a8, const Real &a10, const Real &a12) {
        return (((Real) (46994831769600 * a6 * a6) + (Real) (5482730373120 * a4 * a4) + (Real) (856676620800 * a2 * a2) + (Real) (59227947977932800 * a12 * a12) - 0.106234700e9 * (Real) (a12 * a12) * pow(p, 0.28e2) * pow(x, 0.28e2) - 0.290990700e9 * ((Real) (a8 * a12) + a10 * a10 / 0.2e1 - 0.4e1 * a10 * (Real) a12) * pow(x, 0.24e2) * pow(p, 0.24e2) - 0.347677200e9 * (-0.2e1 * a10 * a10 + (Real) a8 * a10 + (Real) ((a6 - 4 * a8) * a12)) * pow(x, 0.22e2) * pow(p, 0.22e2) - 0.422702280e9 * pow(x, 0.20e2) * ((Real) (a6 - 4 * a8) * a10 + (Real) ((a4 - 4 * a6) * a12) + (Real) (a8 * a8) / 0.2e1) * pow(p, 0.20e2) - 0.524924400e9 * ((Real) (a4 - 4 * a6) * a10 + (Real) ((a2 - 4 * a4) * a12) + (Real) (a8 * (a6 - 2 * a8))) * pow(x, 0.18e2) * pow(p, 0.18e2) - 0.669278610e9 * ((Real) (a2 - 4 * a4) * a10 + (Real) ((a0 - 4 * a2) * a12) + (Real) (a4 * a8) + (Real) (a6 * (a6 - 8 * a8)) / 0.2e1) * pow(x, 0.16e2) * pow(p, 0.16e2) + 0.5103923547340800e16 * a10 * a10 + (Real) (467859658506240 * a8 * a8) - 0.247118256e9 * pow(x, 0.26e2) * (Real) a12 * (a10 - (Real) (2 * a12)) * pow(p, 0.26e2) - 0.1216870200e10 * (-0.4e1 * (Real) a0 * a10 + (Real) (a0 * a8) + (Real) ((a6 - 4 * a8) * a2) + (Real) (a4 * (a4 - 8 * a6)) / 0.2e1) * pow(x, 0.12e2) * pow(p, 0.12e2) - 0.1784742960e10 * (Real) ((a6 - 4 * a8) * a0 + (a4 - 4 * a6) * a2 - 2 * a4 * a4) * pow(x, 0.10e2) * pow(p, 0.10e2) - 0.2868336900e10 * ((Real) ((a4 - 4 * a6) * a0) + (Real) (a2 * (a2 - 8 * a4)) / 0.2e1) * pow(x, 0.8e1) * pow(p, 0.8e1) - 0.5354228880e10 * (Real) ((a2 - 4 * a4) * a0 - 2 * a2 * a2) * pow(x, 0.6e1) * pow(p, 0.6e1) + 0.160626866400e12 * (Real) (a0 * a0) * p * p * x * x + (Real) (321253732800 * a0 * a0) + (Real) (292412286566400 * a6 * a8) - 0.6692786100e10 * pow(x, 0.4e1) * (Real) a0 * (Real) (a0 - 8 * a2) * pow(p, 0.4e1) - 0.882565200e9 * ((Real) (a0 - 4 * a2) * a10 - (Real) (4 * a0 * a12) + (Real) (a2 * a8) + (Real) ((a6 - 4 * a8) * a4) - (Real) (2 * a6 * a6)) * pow(x, 0.14e2) * pow(p, 0.14e2) + (Real) (31329887846400 * a0 + 93989663539200 * a2 + 292412286566400 * a4 + 935719317012480 * a6 + 3062354128404480 * a8 + 34549636320460800 * a12) * a10 + (Real) ((93989663539200 * a0 + 292412286566400 * a2 + 935719317012480 * a4 + 3062354128404480 * a6 + 10207847094681600 * a8) * a12) + (Real) ((856676620800 * a2 + 1713353241600 * a4 + 4112047779840 * a6 + 10965460746240 * a8) * a0) + (Real) ((4112047779840 * a4 + 10965460746240 * a6 + 31329887846400 * a8) * a2) + (Real) ((31329887846400 * a6 + 93989663539200 * a8) * a4)) / (0.1048576e7 / 0.161e3 * a10 * a10 + (0.6144e4 / 0.143e3 * (Real) a0 + 0.8192e4 / 0.65e2 * (Real) a2 + 0.32768e5 / 0.85e2 * (Real) a4 + 0.393216e6 / 0.323e3 * (Real) a6 + 0.524288e6 / 0.133e3 * (Real) a8 + 0.25165824e8 / 0.575e3 * (Real) a12) * a10 + 0.16777216e8 / 0.225e3 * (Real) a12 * (Real) a12 + (0.8192e4 / 0.65e2 * (Real) a0 + 0.32768e5 / 0.85e2 * (Real) a2 + 0.393216e6 / 0.323e3 * (Real) a4 + 0.524288e6 / 0.133e3 * (Real) a6 + 0.2097152e7 / 0.161e3 * (Real) a8) * (Real) a12 + (Real) (a0 * a0) + (0.8e1 / 0.5e1 * (Real) a2 + 0.96e2 / 0.35e2 * (Real) a4 + 0.128e3 / 0.21e2 * (Real) a6 + 0.512e3 / 0.33e2 * (Real) a8) * (Real) a0 + 0.48e2 / 0.35e2 * (Real) a2 * (Real) a2 + (0.128e3 / 0.21e2 * (Real) a4 + 0.512e3 / 0.33e2 * (Real) a6 + 0.6144e4 / 0.143e3 * (Real) a8) * (Real) a2 + 0.256e3 / 0.33e2 * (Real) a4 * (Real) a4 + (0.6144e4 / 0.143e3 * (Real) a6 + 0.8192e4 / 0.65e2 * (Real) a8) * (Real) a4 + 0.4096e4 / 0.65e2 * (Real) a6 * (Real) a6 + 0.196608e6 / 0.323e3 * (Real) a8 * (Real) a8 + 0.32768e5 / 0.85e2 * (Real) a6 * (Real) a8) / p / 0.428338310400e12);
    }

    // same as function presented above, but receives vector of coefficients.
    template<typename Real>
    Real compute_symmetric_fake_log_2cosh(Real x, Real p, const std::vector<Real> &poly_coeffs) {
        Real result = 0.0;

        // Building default values. To create the data structure needed for automatic differentiation, the easiest way is
        // to multiply the default values of coefficients by 1, already containing the derivatives data structure.
        Real one = p/p;
        Real da0 = 0.01 * one;
        Real da2 = 0.0001 * one;
        Real da4 = 0.0 * one;
        Real da6 = 0.0 * one;
        Real da8 = 0.0 * one;
        Real da10 = 0.0 * one;
        Real da12 = 0.0 * one;

        // If x is outside the blending range, we assume the value of the absolute function
        if (x < -2.0/p)
            return -x;
        if (x > 2.0/p)
            return x;

        // Because, aparently, the version I have of Eigen3 does not take care of operations where one of the variables is
        // not constant and not of the automatic diff type, and, also, because I extracted the functions directly from maple,
        // I thought it would be clearer (mainly in numbering) to keep parameters as individual variables instead of an array.
        switch (poly_coeffs.size()) {
            case 0:
                result = compute_symmetric_fake_log_2cosh(x, p, da0, da2, da4, da6, da8, da10, da12);
                break;
            case 1:
                result = compute_symmetric_fake_log_2cosh(x, p, poly_coeffs[0], da2, da4, da6, da8, da10, da12);
                break;
            case 2:
                result = compute_symmetric_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], da4, da6, da8, da10, da12);
                break;
            case 3:
                result = compute_symmetric_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2], da6, da8, da10, da12);
                break;
            case 4:
                result = compute_symmetric_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2], poly_coeffs[3], da8, da10, da12);
                break;
            case 5:
                result = compute_symmetric_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2], poly_coeffs[3], poly_coeffs[4], da10, da12);
                break;
            case 6:
                result = compute_symmetric_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2], poly_coeffs[3], poly_coeffs[4], poly_coeffs[5], da12);
                break;
            case 7:
                result = compute_symmetric_fake_log_2cosh(x, p, poly_coeffs[0], poly_coeffs[1], poly_coeffs[2], poly_coeffs[3], poly_coeffs[4], poly_coeffs[5], poly_coeffs[6]);
                break;
            default:
                throw std::runtime_error("Non recognizable number of extra parameters: " + std::to_string(poly_coeffs.size()));
        }

        if (hasInvalidDerivatives(result) || std::isnan(stripAutoDiff(result)) || std::isinf(stripAutoDiff(result))) {
            std::cerr << "  result: "; std::cerr << result << std::endl;

            std::cerr << "Report derivatives: " << std::endl;
            std::cerr << "  result: "  ; reportDerivatives(std::cerr, result); std::cerr << std::endl; std::cerr << std::endl;
        }

        return result;
    }

    // Smooth version of min function based on polynomial symmetric approximation of absolute function shown above.
    // This function deals with only two entries
    template<typename Real>
    Real exp_smin_symmetric_params(Real a, Real b, Real s, const std::vector<Real> &poly_coeffs) {
        if (s < 1e-10) {
            return std::min(a, b);
        }

        Real epsilon = (a-b) / 2.0;
        Real p = 1.0 / s;
        Real result = (a+b) / 2.0 - compute_symmetric_fake_log_2cosh(epsilon, p, poly_coeffs);

        if (hasInvalidDerivatives(result) || std::isnan(stripAutoDiff(result)) || std::isinf(stripAutoDiff(result))) {
            std::cerr << "  result: "; std::cerr << result << std::endl;
            std::cerr << "  s: "; std::cerr << s << std::endl;
            std::cerr << "  p: "; std::cerr << p << std::endl;
            std::cerr << "  epsilon: "; std::cerr << epsilon << std::endl;

            std::cerr << "Report derivatives: " << std::endl;
            std::cerr << "  result: "  ; reportDerivatives(std::cerr, result); std::cerr << std::endl; std::cerr << std::endl;
            std::cerr << "  s: "  ; reportDerivatives(std::cerr, s); std::cerr << std::endl; std::cerr << std::endl;
            std::cerr << "  p: "  ; reportDerivatives(std::cerr, p); std::cerr << std::endl; std::cerr << std::endl;
            std::cerr << "  epsilon: "  ; reportDerivatives(std::cerr, epsilon); std::cerr << std::endl; std::cerr << std::endl;
        }

        return result;
    }

    // Recursive smooth version of min function based on polynomial symmetric approximation of absolute function shown above.
    // This function deals with an array of values
    template<typename Real>
    Real exp_smin_symmetric_params(const std::vector<Real> &values, Real s, const std::vector<Real> &poly_coeffs) {
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
            result = exp_smin_symmetric_params<Real>(values[0], values[1], s, poly_coeffs);
            return result;
        }

        // Copy values to array to be used in recursive call.
        // At each iteration i, we want to compute recursively, removing the ith from the values array.
        std::vector<Real> recursionValues(values);
        Real lastElement = values[n-1];
        recursionValues.resize(n-1);
        for (size_t i = 0; i < n; i++) {
            Real currentElement = values[i];

            if (i < (n-1)) {
                // exchange current position by last value
                recursionValues[i] = lastElement;
            }

            Real recursiveS = s/10.0;
            Real recursiveMin = exp_smin_symmetric_params<Real>(recursionValues, recursiveS, poly_coeffs);
            result += exp_smin_symmetric_params<Real>(currentElement, recursiveMin, s, poly_coeffs);

            if (i < (n-1)) {
                // put back original element
                recursionValues[i] = currentElement;
            }

            // Check if everyhting is alright!
            if (hasInvalidDerivatives(result) || std::isnan(stripAutoDiff(result)) || std::isinf(stripAutoDiff(result))) {
                std::cerr << "  result: "; std::cerr << result << std::endl;
                std::cerr << "  s: "; std::cerr << s << std::endl;
                std::cerr << "  recursiveMin: "; std::cerr << recursiveMin << std::endl;

                std::cerr << "Report derivatives: " << std::endl;
                std::cerr << "  result: "  ; reportDerivatives(std::cerr, result); std::cerr << std::endl; std::cerr << std::endl;
                std::cerr << "  s: "  ; reportDerivatives(std::cerr, s); std::cerr << std::endl; std::cerr << std::endl;
                std::cerr << "  recursiveMin: "  ; reportDerivatives(std::cerr, recursiveMin); std::cerr << std::endl; std::cerr << std::endl;

                throw std::runtime_error("Error while computing smooth min with a polynomial convex function");
            }
        }
        result /= (Real)n;

        return result;
    }

}

#endif //POLYNOMIALCONVEXBLENDING_H

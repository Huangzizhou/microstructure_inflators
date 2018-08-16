//
// Created by Davi Colli Tozoni on 4/17/18.
//
namespace PiecewiseBlending {
// Functions t, r, h, g, q, fl, flp, f are presented in script doc/KS_piecewise.mw and are explained in document
// doc/blending.pdf
// The general idea is that each coefficient of function f() chooses the curvature of each segment of the whole
// "blending" region. Then, I manually integrated twice the curvatures and took care to guarantee continuity.
    template<typename Real>
    Real t(size_t N, const Real p, const std::vector<Real> &a) {
        size_t PIECES = a.size();

        return (-2 / p * (PIECES - N) / PIECES);
    }

    template<typename Real>
    Real r(size_t N, Real x, const Real p, const std::vector<Real> &a) {
        if (N == 0) {
            return 0;
        } else {
            return a[N - 1] * x - a[N - 1] * t(N - 1, p, a) + r(N - 1, t(N - 1, p, a), p, a);
        }
    }

    template<typename Real>
    Real h(size_t N, Real x, const Real p, const std::vector<Real> &a) {
        if (N == 0) {
            return 0;
        } else {
            return r(N - 1, t(N - 1, p, a), p, a) * (x - t(N - 1, p, a)) +
                   a[N - 1] * pow(t(N - 1, p, a), 0.2e1) / 0.2e1 + a[N - 1] * x * x / 0.2e1 -
                   a[N - 1] * t(N - 1, p, a) * x + h(N - 1, t(N - 1, p, a), p, a);
        }
    }

    template<typename Real>
    Real g(size_t N, Real x, const Real p, const std::vector<Real> &a) {
        size_t PIECES = a.size();

        if (N == 0) {
            return -1;
        } else {
            return h(N, x) / h(PIECES, 0, p, a) - 1;
        }
    }

    template<typename Real>
    Real q(size_t N, Real x, const Real p, const std::vector<Real> &a) {
        size_t PIECES = a.size();

        if (N == 0) {
            return -x + t(0, p, a);
        } else {
            Real zero = 0.0 * p;
            Real multiplier = 0.1e1 / h(PIECES, zero, p, a);
            Real firstTerm = (pow(x, 0.3e1) - pow(t(N - 1, p, a), 0.3e1)) * a[N - 1] / 0.6e1;
            Real secondTermMultiplier = (x * x - pow(t(N - 1, p, a), 0.2e1));
            Real secondTermLhs = r(N - 1, t(N - 1, p, a), p, a) / 0.2e1;
            Real secondTermRhs = a[N - 1] * t(N - 1, p, a) / 0.2e1;
            Real secondTerm = secondTermMultiplier * (secondTermLhs - secondTermRhs);
            Real thirdTerm = (h(N - 1, t(N - 1, p, a), p, a) - r(N - 1, t(N - 1, p, a), p, a) * t(N - 1, p, a) +
                              a[N - 1] * pow(t(N - 1, p, a), 0.2e1) / 0.2e1);
            Real rhs = -x + t(N - 1, p, a) + q(N - 1, t(N - 1, p, a), p, a);
            return multiplier * (firstTerm + secondTerm + (x - t(N - 1, p, a)) * thirdTerm) + rhs;
        }
    }

    template<typename Real>
    Real fl(size_t N, Real x, const Real p, const std::vector<Real> &a) {
        if (N == 0) {
            return -x;
        } else {
            return q(N, x, p, a) - t(0, p, a) - q(0, t(0, p, a), p, a);
        }
    }

    template<typename Real>
    Real flp(Real x, const Real p, const std::vector<Real> &a) {
        Real result = 0.0;
        size_t PIECES = a.size();

        if (x <= -2 / p) {
            result = fl(0, x, p, a);
        } else if (x <= -2 / p * (PIECES - 1) / PIECES) {
            result = fl(1, x, p, a);
        } else if (x <= -2 / p * (PIECES - 2) / PIECES) {
            result = fl(2, x, p, a);
        } else if (x <= -2 / p * (PIECES - 3) / PIECES) {
            result = fl(3, x, p, a);
        } else if (x <= -2 / p * (PIECES - 4) / PIECES) {
            result = fl(4, x, p, a);
        } else if (x <= -2 / p * (PIECES - 5) / PIECES) {
            result = fl(5, x, p, a);
        } else if (x <= -2 / p * (PIECES - 6) / PIECES) {
            result = fl(6, x, p, a);
        } else if (x <= -2 / p * (PIECES - 7) / PIECES) {
            result = fl(7, x, p, a);
        } else if (x <= -2 / p * (PIECES - 8) / PIECES) {
            result = fl(8, x, p, a);
        }

        return result;
    }

    template<typename Real>
    Real f(Real x, const Real p, const std::vector<Real> &a) {
        Real result;
        if (x <= 0) {
            result = flp(x, p, a);
        } else {
            Real negativeX = -1.0 * x;
            result = flp(negativeX, p, a);
        }

        return result;
    }

}

namespace SD {

    // This function corresponds to an approximation of the absolute function. The name is refers to the fact that this is
    // a piecewise symmetric function. Each coefficient corresponds to the curvature expected in each segment of the blending
    // region. A more detailed explanation is presented doc/blending.pdf).
    // Function code was based on code extracted from Maple script doc/KS_piecewise.mw
    template<typename Real>
    Real compute_piecewise_fake_log_2cosh(Real x, Real p, const std::vector<Real> poly_coeffs) {
        Real result = 0.0;

        // Building default values
        Real one = p / p;
        Real da0 = 0.01 * one;

        // If x is outside the blending range, we assume the value of the absolute function
        if (x < -2.0 / p)
            return -x;
        if (x > 2.0 / p)
            return x;

        // Call function f(), responsible for the computation of the approximate absolute value
        std::vector<Real> local_coeffs;
        switch (poly_coeffs.size()) {
            case 0:
                local_coeffs.push_back(da0);
                result = PiecewiseBlending::f(x, p, local_coeffs);
                break;
            default:
                if (poly_coeffs.size() > 8)
                    throw std::runtime_error("Non recognizable number of extra parameters: " + std::to_string(poly_coeffs.size()));

                result = PiecewiseBlending::f(x, p, poly_coeffs);
                break;
        }

        return result;
    }

    // Smooth version of min function based on piecewise approximation of absolute function shown above.
    // This function deals with only two entries
    template<typename Real>
    Real exp_smin_piecewise_params(Real a, Real b, Real s, const std::vector<Real> &poly_coeffs) {
        if (s < 1e-20) {
            return std::min(a, b);
        }

        Real epsilon = (a - b) / 2.0;
        Real p = 1.0 / s;
        Real result;


        result = (a + b) / 2.0 - compute_piecewise_fake_log_2cosh(epsilon, p, poly_coeffs);

        return result;
    }

     // Recursive smooth version of min function based on piecewise symmetric approximation of absolute function shown above.
    // This function deals with an array of values
    template<typename Real>
    Real exp_smin_piecewise_params(const std::vector<Real> &values, Real s, const std::vector<Real> &poly_coeffs) {
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
            result = SD::exp_smin_piecewise_params<Real>(values[0], values[1], s, poly_coeffs);
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
            Real recursiveMin = exp_smin_piecewise_params<Real>(recursionValues, recursiveS, poly_coeffs);
            result += exp_smin_piecewise_params<Real>(currentElement, recursiveMin, s, poly_coeffs);

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

                throw std::runtime_error("Error while computing smooth min with a piecewise function");
            }
        }
        result /= (Real) n;

        return result;
    }
}


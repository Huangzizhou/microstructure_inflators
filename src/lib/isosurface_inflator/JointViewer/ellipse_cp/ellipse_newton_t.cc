////////////////////////////////////////////////////////////////////////////////
// ellipse_newton.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Solves for the closest point on an axis-aligned ellipse with major axis
//  a and minor axis b <= a:
//  e(t) = [a cos(t), b sin(t)]
//  (Note: the major axis is aligned with the first coordinate).
//  Applies Newton's method with a carefully chosen starting point that
//  guarantees convergence.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/31/2016 14:48:46
////////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <fstream>

#include <iostream>
#include <tuple>
#include <ratio>
#include <cassert>
#include <stdexcept>
#include <cmath>

#ifndef OUTPUT_FILE
#define OUTPUT_FILE 1
#endif

// For ellipse with major axis 1, minor axis b:
//  e(t) = [cos(t), b sin(t)]
// and query point in positive octant
// Returns closest point's ellipse parameter, t in [0, pi/2], and distance:
// (t, dist)
template<typename Real, bool verbose = false, typename TOL = std::ratio<1,long(1e9)>>
std::tuple<Real,Real,size_t,size_t>
closestPoint(Real b, Real x, Real y) {
    if (!((x >= 0) && (y >= 0)))
        throw std::runtime_error("Query point in incorrect quadrant.");

    // Number of Newton iterations performed
    size_t it = 0;
    // Number of backtracking iterations
    size_t backtrack_it = 0;

    Real tol = Real(TOL::num)/ TOL::den;

    // Evaluating close to the center of the ellipse ((x, y) approx 0) will lead
    // to numerical issues. The center, and more generally all points near the
    // minor axis, will have a closest point at the same x value.
    if (x < 1e-8) { // ~= sqrt(epsilon)
        // ellipse y coordinate at x:
        //      ey = b sqrt(1 - x * x) ~= b to double precision
        Real ey = b;
        return {acos(x), y - ey, it, backtrack_it};
    }

    // Determine if we're inside or outside the ellipse
    // When f < 0, the point is inside and f > 0 the point is outside:
    Real xSq = x * x, ySq = y * y, bSq = b * b;
    Real f = xSq + ySq / bSq - 1;

    // There is a singularity inside the ellipse, lying on the major axis that
    // causes slow convergence in its vicinity. However, we have an explicit
    // formula for the closest point on the major axis:
    if ((f < 0) && (y < b * 1e-6)) {
        if (x < (1 - bSq)) {
            // Closest point is on ellipse interior
            Real t = acos(x / (1 - bSq)); // Argument guaranteed < 1
            Real distSq = (bSq * (xSq + (bSq - 1))) / (bSq - 1);
            return {t, copysign(sqrt(distSq), f), it, backtrack_it};
        }
        else {
            // Closest point is at ellipse endpoint t = 0
            return {0, x - 1, it, backtrack_it};
        }
    }

    struct TemporaryCache {
        TemporaryCache() { }
        TemporaryCache(Real b, Real x, Real y, Real t)         { updateIntermediates(b, x, y, t); }
        TemporaryCache(Real b, Real x, Real y, Real c, Real s) { updateIntermediates(b, x, y, c, s); }

        // Cached values for intermediates needed in all computations
        Real cos_t, sin_t, bsin_t;
        Real xdist, ydist;

        // Cached value for derivatives (so that values computed in line search can
        // be reused in next Newton iteration.)
        Real dsqdist, d2sqdist;

        // From angle
        void updateIntermediates(Real b, Real x, Real y, Real t) { updateIntermediates(b, x, y, cos(t), sin(t)); }

        // From cos, sin
        void updateIntermediates(Real b, Real x, Real y, Real c, Real s) {
            cos_t = c, sin_t = s;
            bsin_t = b * sin_t;
            xdist = x - cos_t;
            ydist = y - bsin_t;
        }

        Real evalSqDist() const { return xdist * xdist + ydist * ydist; };

        // ||[x,y] - e(t)||^2 and its first and second derivatives.
        // Note: computes 0.5 * derivatives (Newton step is invariant to scaling)
        // Must be called after updateIntermediates!
        void updateFirstDerivative(Real b)                    {  dsqdist = xdist * sin_t - ydist * b * cos_t; }
        void updateSecondDerivative(Real bSq, Real x, Real y) { d2sqdist = (bSq - 1) * (cos_t * cos_t - sin_t * sin_t) + x * cos_t + y * bsin_t; }
    };

    // In the general case, we run Newton's method from an initial value
    // carefully chosen to lie within an interval that contains the optimum and
    // over which the distance function is convex. Then our line search ensures
    // that we remain in this interval, guaranteeing convergence.
    //
    // It turns out that all points counterclockwise (with greater t value) from
    // the closest ellipse point satisfy this criterion. We can pick such a
    // starting point as follows:
    //      t_0 = atan(y/(bx)) if (y, x) outside the ellipse
    //            acos(x)      if (y, x) inside  the ellipse
    //
    // Both of these initial parameters converge to the optimum as (x, y) moves
    // toward the ellipse, so they are reasonable initial guesses.
#if 1
    // Try the atan guess:
    // atan2(y, b x) ==> Normalize[{b * x, y}]
    Real invNorm = 1 / sqrt(ySq + bSq * xSq);
    TemporaryCache c_atan(b, x, y, b * x * invNorm, y * invNorm);

    // Try the acos guess:
    //  acos(x) ==> {x, sqrt(1 - x^2)}
    TemporaryCache c_acos;
    if (x < 1) c_acos.updateIntermediates(b, x, y, x, sqrt(1 - xSq));
    else       c_acos.updateIntermediates(b, x, y, 1,             0);

    // At least one of these guesses is guaranteed to be in the convex interval.
    // Pick the closest one if it is in the interval (generally the case), and
    // resort to the other if not.
    // (Never computes more second derivatives/trig functions than necessary)
    TemporaryCache c;
    Real t;
    if (c_atan.evalSqDist() < c_acos.evalSqDist()) {
        c_atan.updateSecondDerivative(bSq, x, y);
        if (c_atan.d2sqdist <= 0) {
            c_acos.updateSecondDerivative(bSq, x, y);
            c = c_acos;
            t = (x >= 1) ? 0.0 : acos(x);
        }
        else {
            c = c_atan;
            t = atan2(y, b * x);
        }
    }
    else {
        c_acos.updateSecondDerivative(bSq, x, y);
        if (c_acos.d2sqdist <= 0) {
            c_atan.updateSecondDerivative(bSq, x, y);
            c = c_atan;
            t = atan2(y, b * x);
        }
        else {
            c = c_acos;
            t = (x >= 1.0) ? 0.0 : acos(x);
        }
    }

    assert(c.d2sqdist > 0);

    // Note: the midpoint of these two guesses is often a better choice but
    // requires more trig functions to compute, which undoes the performance
    // gains in practice.
#else
    Real t = f > 0 ? atan2(y, b * x) : acos(std::min<Real>(x, 1.0));
    TemporaryCache c(b, x, y, t);
    c.updateSecondDerivative(bSq, x, y);
#endif
    // c already has second derivative information for t (but not first)
    c.updateFirstDerivative(b);

    // Newton iteration
    // Loop invariant:
    // Upon entry, c.dsqdist and c.d2sqdist hold derivatives for current t.
    Real delta_t = 0;
    while (true) { // Convergence check after delta_t is computed...
        if (verbose) std::cerr << t << '\t' << sqrt(c.evalSqDist()) << '\t' << c.dsqdist << '\t' << c.d2sqdist << std::endl;
        ++it;
        delta_t = -c.dsqdist / c.d2sqdist;
        if (std::abs(delta_t) < tol) break;
        if (c.d2sqdist <= 0) {
            std::cerr << "x: " << x << ", y: " << y << std::endl;
            std::cerr << t << '\t' << sqrt(c.evalSqDist())
                      << '\t' << c.dsqdist
                      << '\t' << c.d2sqdist
                      << std::endl;
            assert(false);
        }

        // Clamp overly large steps to the correct quadrant
        // This case happens, e.g., for points near the center of a circle.
        delta_t = std::max<Real>(delta_t, -t);
        delta_t = std::min<Real>(delta_t, M_PI / 2 - t);

        // Converged?
        if (std::abs(delta_t) <= tol) {
            t += delta_t; // accept the (tiny) step
            c.updateIntermediates(b, x, y, t);
            break;
        }

        // Backtracking line search to stay in convex interval
        Real alpha = 1.0;
        Real t_next;
        while (true) {
            t_next = t + alpha * delta_t;
            c.updateIntermediates(b, x, y, t_next); // cost: 1 sin,cos eval
            c.updateSecondDerivative(bSq, x, y);

            if (c.d2sqdist >= 0) break; // still in convex region.

            // Overstepped into concave region; backtrack
            alpha *= 0.75;
            ++backtrack_it;
        }

        // second derivative is already current, but first must be updated
        c.updateFirstDerivative(b);
        t = t_next;
    }; // (convergence criterion based on non-backtracked Newton step size)

    return {t, copysign(sqrt(c.evalSqDist()), f), it, backtrack_it};
}

// For ellipse with major axis a; reduce to major axis 1 case.
//  e(t) = [a cos(t), b sin(t)]
// and query point in positive octant
// Returns closest point's ellipse parameter, t in [0, pi/2]
template<typename Real, bool verbose = false, typename TOL = std::ratio<1,long(1e9)>>
std::tuple<Real,Real, size_t, size_t> 
closestPoint(Real a, Real b, Real x, Real y) {
    if ((a < 0) || (b < 0) || (a < b))
        throw std::runtime_error("Invalid ellipse axis lengths (must be positive and must satisfy a < b).");
    auto result = closestPoint<Real, verbose, TOL>(b/a, x/a, y/a);
    std::get<1>(result) *= a;
    return result;
}

int main(int argc, const char *argv[]) {
    if (!((argc == 3) || (argc == 5))) {
        std::cerr << "usage: ./ellipse_newton a b [x y]" << std::endl;
        exit(-1);
    }

    double a = std::stod(argv[1]),
           b = std::stod(argv[2]);

    double aspect = std::min(a / b, 3.0);
    size_t gridSizeX = 1024,
           gridSizeY = ceil(gridSizeX / aspect);
    if (gridSizeY < 2) throw std::runtime_error("Invalid aspect ratio");
    double gridSpacing = 1.5 * a / (gridSizeX - 1);
    // For now, always sample square domain--makes creating animations easier.
    gridSizeY = gridSizeX;

    std::cout << std::setprecision(19);
    std::cerr << std::setprecision(19);

#if OUTPUT_FILE
    std::ofstream iterateFile("iterates.txt");
    iterateFile << std::setprecision(19);
#endif

    if (argc > 3) {
        double x = std::stod(argv[3]),
               y = std::stod(argv[4]);
        auto cp = closestPoint<double, true>(a, b, x, y);
        std::cout << std::get<0>(cp) << "\t" << std::get<1>(cp) << std::endl;
    }
    else {
        for (size_t xi = 0; xi < gridSizeX; ++xi) {
            for (size_t yi = 0; yi < gridSizeY; ++yi) {
                double x = xi * gridSpacing,
                       y = yi * gridSpacing;
                auto cp = closestPoint<double>(a, b, x, y);
#if OUTPUT_FILE
                iterateFile << x
                    << '\t' << y
                    << '\t' << std::get<0>(cp)
                    << '\t' << std::get<1>(cp)
                    << '\t' << std::get<2>(cp)
                    << '\t' << std::get<3>(cp)
                    << '\n';
#endif
            }
#if OUTPUT_FILE
            iterateFile  << '\n'; // gnuplot wants blank lines between x values
#endif
        }
    }
    return 0;
}

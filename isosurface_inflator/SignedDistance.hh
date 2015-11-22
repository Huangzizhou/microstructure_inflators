////////////////////////////////////////////////////////////////////////////////
// SignedDistance.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Signed distance function utilities. These functions are templated to
//      support auto-differentiation via adept. 
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/23/2015 14:59:39
////////////////////////////////////////////////////////////////////////////////
#ifndef SIGNEDDISTANCE_HH
#define SIGNEDDISTANCE_HH

#include "InflatorTypes.hh"
#include <vector>

namespace SD {

template<typename Real>
Real clamp(Real x, Real a, Real b) { return std::min<Real>(std::max<Real>(x, a), b); }
template<typename Real>
Real mix(Real x, Real y, Real a) { return x * (1 - a) + y * a; }

// Min of a collection of values.
template<typename Real>
Real min(const std::vector<Real> &values) {
    auto min_it = std::min_element(values.begin(), values.end());
    if (min_it == values.end()) return 1e5;
    return *min_it;
}

// exponential smooth min of two values (k = 32);
template<typename Real>
Real exp_smin(Real a, Real b, Real k = 32)
{
    Real res = exp(-k * a) + exp(-k * b);
    return -log(res) / k;
}

// exponential smooth min one or more values (k = 32);
template<typename Real>
Real exp_smin(const std::vector<Real> &values, Real k = 32)
{
    Real res = 0;
    for (Real v : values) {
        res += exp(-k * v);
    }
    return -log(res) / k;
}

// Minimum max(a, b) - min(a, b) such that
// min(a, b) - exp_smin(a, b) < tol
// (Note: exp_smin always under-estimates, so this is a bound on the absolute
// error of the min approximation).
template<typename Real>
Real exp_smin_radius(Real k = 20, Real tol = 1e-3) {
    return -log(exp(k * tol) - 1.0) / k;
}

template<typename Real>
Real poly_smin(Real a, Real b, Real k = 0.2) {
    Real h = clamp(0.5 + 0.5 * (b - a) / k, Real(0.0), Real(1.0));
        return mix(b, a, h) - k * h * (1.0 - h);
}

// power smooth min (k = 8);
template<typename Real>
Real pow_smin( Real a, Real b, Real k = 5 )
{
    a = pow( a, k ); b = pow( b, k );
    return pow( (a*b)/(a+b), 1.0/k );
}

template<typename Real>
Real smax(Real a, Real b, Real k = 0.2) {
    return -poly_smin(-a, -b, k);
}

////////////////////////////////////////////////////////////////////////////////
// Primitives
////////////////////////////////////////////////////////////////////////////////
namespace Primitives {

template<typename Real>
struct Sphere {
    Sphere() { }
    Sphere(const Point3<Real> &c, Real r) { set(c, r); }
    void set(const Point3<Real> &c, Real r) { m_c = c; m_r = r; }

    Real signedDistance(const Point3<Real> &p) const {
        return sqrt((p - m_c).squaredNorm()) - m_r;
    }
private:
    Point3<Real> m_c;
    Real m_r;
};

template<typename Real>
struct Box {
    Box() { }
    Box(const BBox<Point3<Real>> &box) : m_box(box) { }
    void set(const BBox<Point3<Real>> &box) { m_box = box; }
    Real signedDistance(const Point3<Real> &p) const {
        // Box is union of 6 halfspaces
        Vector3<Real> d = (m_box.minCorner - p).cwiseMax(p - m_box.maxCorner);
        return std::min(d.maxCoeff(), Real(0)) + d.cwiseMax(Vector3<Real>::Zero()).norm();
    }
private:
    BBox<Point3<Real>> m_box;
};

template<typename Real>
struct ConicalFrustum {
    ConicalFrustum() { }
    ConicalFrustum(const Point3<Real> &a, const Point3<Real> &b, Real ra, Real rb) { set(a, b, ra, rb); }

    void set(const Point3<Real> &a, const Point3<Real> &b, Real ra, Real rb) {
        m_a = a;
        m_ra = ra; m_rb = rb;
        m_axis = b - a;
        m_axisLenSq = m_axis.squaredNorm();
        m_axisLen = sqrt(m_axisLenSq);
    }

    Real signedDistance(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_a);
        Real frac_along = v.dot(m_axis) / m_axisLenSq;
        Real dist_axial = m_axisLen * (fabs(frac_along - 0.5) - 0.5);

        Vector3<Real> v_perp = v - frac_along * m_axis;
        Real r = m_ra + clamp(frac_along, Real(0.0), Real(1.0))  * (m_rb - m_ra);
        Real dist_perp = sqrt(v_perp.squaredNorm()) - r;

        Vector2<Real> posDists(std::max(dist_perp, Real(0.0)), std::max(dist_axial, Real(0.0)));
        return sqrt(posDists.squaredNorm()) + std::min(std::max(dist_perp, dist_axial), Real(0.0));
    }

private:
    Point3<Real> m_a;
    Vector3<Real> m_axis;
    Real m_axisLen, m_axisLenSq, m_ra, m_rb;
};

template<typename Real>
struct Cylinder {
    Cylinder() { }
    Cylinder(const Point3<Real> &a, const Point3<Real> &b, Real r) { set(a, b, r); }

    void set(const Point3<Real> &a, const Point3<Real> &b, Real r) {
        m_a = a;
        m_axis = b - a;
        m_r = r;
        m_axisLenSq = m_axis.squaredNorm();
        m_axisLen = sqrt(m_axisLenSq);
    }

    Real signedDistance(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_a);
        Real frac_along = v.dot(m_axis) / m_axisLenSq;
        Real dist_axial = m_axisLen * (std::abs(frac_along - 0.5) - 0.5);

        Vector3<Real> v_perp = v - frac_along * m_axis;
        Real dist_perp = sqrt(v_perp.squaredNorm()) - m_r;

        Vector2<Real> posDists(std::max(dist_perp, 0.0), std::max(dist_axial, 0.0));
        return sqrt(posDists.squaredNorm()) + std::min(std::max(dist_perp, dist_axial), 0.0);
    }

private:
    Point3<Real> m_a;
    Vector3<Real> m_axis;
    Real m_axisLen, m_axisLenSq, m_r;
};

// InflatedEdge "Primitive": edge geometry for an inflated wire mesh.
// The edge geometry consists of a "cylinder" with sphere endpoints. To support
// per-vertex thickness, the spheres can be of different radii, in which case
// the cylinder part is really a conical frustum (linearly interpolated
// radius). The conical frustum is created so that it joins smoothly with the
// sphere endpoints (tangent at the intersection--continuous normal).
template<typename Real>
class InflatedEdge {
public:
    InflatedEdge() { }
    InflatedEdge(const Point3<Real> &p1, const Point3<Real> &p2, Real r1, Real r2)
    { set(p1, p2, r1, r2); }

    void set(const Point3<Real> &p1, const Point3<Real> &p2, const Real r1, const Real r2) {
        m_c1 = p1, m_c2 = p2;
        m_r1 = r1, m_r2 = r2;

        m_axisUnit = p2 - p1;
        m_axisLength = sqrt(m_axisUnit.squaredNorm());
        m_axisUnit /= m_axisLength;

        m_sinTheta = -(r2 - r1) / m_axisLength;
        m_theta = asin(m_sinTheta);
        m_cosTheta = sqrt(1 - m_sinTheta * m_sinTheta);

        m_edgeLength = m_axisLength * m_cosTheta;
    }

    // Additional real type to support automatic differentiation wrt. p
    // even when the class's real type doesn't support autodiff.
    template<typename Real2>
    Real2 signedDistance(const Point3<Real2> &p) const {
        // Note: we must cast the vector-type member variables to Real2 types
        // for the automatic differentation case. There's a chance this will
        // make the Real2 == Real case less efficient; we can add some
        // SFINAE vodoo to optimize this if it's an issue.

        Vector3<Real2> v(p - m_c1.template cast<Real2>());
        Real2 v_parallelComponent = v.dot(m_axisUnit.template cast<Real2>());
        Vector3<Real2> v_perp = v - v_parallelComponent * m_axisUnit.template cast<Real2>();
        Real2 v_perpComponent = sqrt(v_perp.squaredNorm());
        // Rotate so that conical frustum surface is horizontal
        Real2 x = m_cosTheta * v_parallelComponent - m_sinTheta * v_perpComponent;

        // Closest surface is sphere 1
        if (x < 0) return sqrt(v.squaredNorm()) - m_r1;

        // Closest surface is the conical frustum part (the closest edge of
        // which is horizontal in rotated (x, y) coordinates).
        Real2 y = m_sinTheta * v_parallelComponent + m_cosTheta * v_perpComponent;
        if (x < m_edgeLength) return y - m_r1;

        // Closest surface is sphere 2 
        return sqrt((p - m_c2.template cast<Real2>()).squaredNorm()) - m_r2;
    }

    Point3<Real> closestMedialAxisPoint(const Point3<Real> &p) const {
        Real parallelDist = (p - m_c1).dot(m_axisUnit);
        // Clamp 'parallelDist' to stay in the line segment (0, m_axisLength)
        return m_c1 + clamp(parallelDist, Real(0), m_axisLength) * m_axisUnit;
    }

    // Signs are based on whether the overlap with the enpoint sphere increases
    // (positive) or decreases (negative)
    Real angleAtP1() const { return  m_theta; }
    Real angleAtP2() const { return -m_theta; }

private:
    Point3<Real> m_c1, m_c2;
    Real m_r1, m_r2;

    Real m_theta, m_sinTheta, m_cosTheta, m_axisLength, m_edgeLength;
    Vector3<Real> m_axisUnit;
};


} // End of namespace: Primitives

}

#endif /* end of include guard: SIGNEDDISTANCE_HH */

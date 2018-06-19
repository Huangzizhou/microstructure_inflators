#include <AntTweakBar.h>
#include <MeshFEM/Future.hh>
#include <MeshFEM/utils.hh>
#include <MeshFEM/colors.hh>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <functional>
#include <sstream>

#if defined(__APPLE__)
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

#include "quaternions.hh"

// #include "../VCGSurfaceMesherImpl.hh"
#include "../IGLSurfaceMesherMC.hh"
#include "../InflatorTypes.hh"
#include "../SignedDistance.hh"

#include "ConvexHull.hh"

#include "ellipse_cp/ellipse_newton.hh"
#include <MeshFEM/filters/gen_grid.hh>
#include <MeshFEM/MSHFieldWriter.hh>

void drawArrow(const Point3D &tail,
               const Point3D &tip, Real radius,
               const RGBColorf &color = RGBColorf(1.0, 0.0, 0.0));

typedef enum {SMOOTH_NONE, SMOOTH_FULL, SMOOTH_ELLIPSECP, SMOOTH_RAYCAST, SMOOTH_HULL} SmoothingMode;

struct GeometryParams {
    float radius1 = 0.05, radius2 = 0.20, radius3 = 0.20;
    float alpha = 0.85;
    float edgeLen1 = 0.75, edgeLen2 = 0.75;

    // float radius1 = 0.20, radius2 = 0.05, radius3 = 0.5;
    // float alpha = 0.8;
    int mcGridSize = 99;

    bool smooth = false;
    SmoothingMode smoothingMode = SMOOTH_NONE;
    float smoothness = 0.01;
    float hullModulationWidth = 0.005;

    float isolevel = 0.0;

    // For determining if we need to remesh only.
    bool operator==(const GeometryParams &b) const {
        return (mcGridSize          == b.mcGridSize) &&
               (radius1             == b.radius1) &&
               (radius2             == b.radius2) &&
               (radius3             == b.radius3) &&
               (alpha               == b.alpha) &&
               (edgeLen1            == b.edgeLen1) &&
               (edgeLen2            == b.edgeLen2) &&
               (smoothingMode       == b.smoothingMode) &&
               (smoothness          == b.smoothness) &&
               (hullModulationWidth == b.hullModulationWidth) &&
               (isolevel            == b.isolevel);
    }
};

// Represents a conic section in (t, z) coordinates
// Intersection curve (conic section) in (t, z) variables
// ((t - c) / a)^2 + (z / b)^2 = 1
//
// z = +/-b sqrt(1 - ((t - c)/a)^2)
//   = +/-(b/a) sqrt(a^2 - (t - c)^2)
// t = c +/- a sqrt(1 - (z/b)^2)
//   = c +/- sqrt(a^2 - (a/b)^2 z^2)
template<typename Real>
struct ConicSection {
public:
    // Choose between ellipse, parabola, and hyperbola based on the
    // magnitude of b/a
    // (which is computed in a more numerically stable way than a or b)
    ConicSection() : m_a(0), m_b_div_a(0), m_c(0) { }
    ConicSection(Real a, Real b_div_a, Real a_plus_c)
                : m_a(a), m_b_div_a(b_div_a), m_a_plus_c(a_plus_c), m_c(a_plus_c - a) {
        assert(std::abs(m_b_div_a) <= 1 + 1e-8);
        m_b_div_a = std::min<Real>(m_b_div_a, 1);
    }

    Real z(Real tval) const {
        tval -= m_c;
        return m_b_div_a * sqrt((m_a + tval) * (m_a - tval));
    }

    Real t(Real zval) const {
        Real sgn = (m_a < 0) ? -1.0 : 0.0;
        return m_c + sgn * sqrt(m_a * m_a - (zval * zval) / (m_b_div_a * m_b_div_a));
    }

    std::pair<Real, Real> pointAtAngle(Real theta) const {
        return { m_a * cos(theta) + m_c,
                 m_a * m_b_div_a * sin(theta) };
    }

    std::pair<Real, Real> normalAtPoint(const std::pair<Real, Real> &tz) const {
        Real t = tz.first, z = tz.second;
        Real grad_t = t - m_c;
        Real grad_z = (z / m_b_div_a) / m_b_div_a;
        Real invNorm = 1.0 / sqrt(grad_t * grad_t + grad_z * grad_z);
        return { grad_t * invNorm, grad_z * invNorm };
    }

    std::pair<Real, Real> closestPtCoords(const std::pair<Real, Real> &tz) const {
        Real t = tz.first, z = tz.second;
        // Ellipse case:
        // Translate to center
        t -= m_c;
        // Reduce to the positive quadrant with reflections
        Real tSign = 1.0, zSign = 1.0;
        if (t < 0) { tSign = -1.0; t *= -1; }
        if (z < 0) { zSign = -1.0; z *= -1; }

        Real cp_t, cp_z, dist;
        size_t dummy;
        std::tie(cp_t, cp_z, dist, dummy, dummy) =
            ellipseClosestPoint(m_b_div_a, t / m_a, z / m_a);

        // Undo scalings, reflections and translations
        cp_t *= m_a * tSign;
        cp_z *= m_a * zSign;
        cp_t += m_c;
        return { cp_t, cp_z };
    }

private:
    Real m_a, m_b_div_a, m_a_plus_c, m_c;
};

// Polynomial whose zero isosurface is a cone/cylinder along the x axis with
// radius "b" at the origin and "slope" m.
template<typename _Real = double>
struct ConeQuadric {
    using Real = _Real;
    using XF = std::function<Point3D(const Point3D &p)>;

    ConeQuadric(Real m, Real b, const XF &xf = nullptr) : m(m), b(b), xform(xf) { }

    // Cone/cylinder coinciding with inflated geometry for **unit length** edge
    // lying in xy-plane, making an angle "alpha" with the x axis and with
    // endsphere radii r1 and r2.
    static ConeQuadric QuadricXY(Real alpha, Real r1, Real r2) {
        Real sin_th = r2 - r1;
        Real cos_th = sqrt(1 - sin_th * sin_th);
        Real x1 = -r1 * sin_th;
        Real l = 1 - r2 * sin_th + r1 * sin_th;
        Real m = (r2 - r1) * cos_th / l;
        Real b = r1 * cos_th;
        // std::cout << "b:\t" << b << std::endl;
        return ConeQuadric(m, b,
            [=](const Point3D &p) { return Point3D(p[0] * cos(alpha) + p[1] * sin(alpha) - x1,
                                                  -p[0] * sin(alpha) + p[1] * cos(alpha), p[2]); }
            );
    }

    static ConeQuadric FromCanonicalParams(Real alpha, Real r2) {
        return QuadricXY(alpha, 1.0, r2);
    }

    // Cone/cylinder coinciding with inflated geometry for edge
    // lying in xy-plane, making an angle "alpha" with the x axis.
    // The edge has length "l" and endsphere radii r1 and r2.
    // (Convert to equivalent unit length edge)
    static ConeQuadric FromXYEdge(Real alpha, Real l, Real r1, Real r2) {
        Real sin_th = (r2 - r1) / l;
        return QuadricXY(alpha, r1, r1 + sin_th);
    }

    Real signedDistance(Point3D p) const {
        if (xform) p = xform(p);
        return p[1] * p[1] + p[2] * p[2] - (m * p[0] + b) * (m * p[0] + b);
    }
    BBox<Point3D> boundingBox() const { return BBox<Point3D>(Point3D(-1, -1, -1), Point3D(1, 1, 1)); }

    Real b, m;
    XF xform;
};

std::pair<ConeQuadric<>, ConeQuadric<>> getQuadrics();
//     * c3
// e2 /
//   / a
//c1+----* c2
//    e1
template<typename Real = double>
struct BinaryJoint {
    using Vec = Vector3<Real>;
    using  Pt =  Point3<Real>;

    BinaryJoint() { }

    BinaryJoint(const Pt &c1, const Pt &c2, const Pt &c3,
                Real r1, Real r2, Real r3) {
        set(c1, c2, c3, r1, r2, r3);
    }

    void set(const Pt &c1, const Pt &c2, const Pt &c3,
             Real r1, Real r2, Real r3) {
        Vec e1 = c2 - c1,
            e2 = c3 - c1;
        Real l1 = e1.norm(),
             l2 = e2.norm();
        e1 /= l1;
        e2 /= l2;

        Real cosAlpha = e1.dot(e2);
        // Unit vector perpendicular to edge pair
        Vec z = e1.cross(e2);
        Real sinAlpha = z.norm();
        if ((sinAlpha < 1e-6) && (cosAlpha > 0)) {
            m_degenerate = true; // Signal a degenerate joint: smoothing will be disabled
            z = { 0, 0, 1 }; // Arbitrary
            sinAlpha = 1;
        }

        // Construct a unit vector (as measured in 3D) pointing in the
        // intersection plane's z direction.
        z /= sinAlpha;

        Real alpha = atan2(sinAlpha, cosAlpha);
        Real c_r2, c_r3;
        std::tie(c_r2, c_r3) =
            m_canonicalParamsFromParams(l1, l2, r1, r2, r3);

        // Construct plane normal with coefficients nx, ny in frame (e1, e1perp)
        Vec e1perp = z.cross(e1);
        Real nx, ny, planeIsolevel;

        Real      r2_2_minus_r2 = (2.0 - c_r2) * c_r2;
        Real sqrt_r2_2_minus_r2 = sqrt(r2_2_minus_r2);
        Real sqrtTerm = sqrt(((2.0 - c_r3) * c_r3)) / sqrt_r2_2_minus_r2;
        nx = cosAlpha - sqrtTerm;
        ny = sinAlpha;
        planeIsolevel = (1.0 - c_r3) - (1.0 - c_r2) * sqrtTerm;
        Real planeNormalMagSq = nx * nx + ny * ny;
        Real planeNormalMag   = sqrt(planeNormalMagSq);
        m_iplaneNormal = (nx * e1 + ny * e1perp) / planeNormalMag;

        // t direction: tangent vector to intersection plane,
        // perp to z (the other tangent vector)
        Vec t = (ny * e1 - nx * e1perp) / planeNormalMag;
        m_iplaneOrigin  = nx * e1 + ny * e1perp;
        m_iplaneOrigin *= r1 * planeIsolevel / planeNormalMagSq;
        m_iplaneOrigin += c1;

        assert(std::abs(t.squaredNorm() - 1.0) < 1e-8);
        assert(std::abs(z.squaredNorm() - 1.0) < 1e-8);

        // Our (t, z) parametrization of the intersection plane uses "canonical"
        // units in which r1 is 1.
        m_dt   = t / r1;
        m_dz   = z / r1;
        m_d_dt = t * r1;
        m_d_dz = z * r1;

        // Intersection curve (conic section) in (t, z) variables
        // ((t - c) / a)^2 + (z / b)^2 = 1
        Real den = planeNormalMagSq * r2_2_minus_r2 - ny * ny;
        Real num = planeIsolevel * (1 - c_r2) - nx;

        Real a = planeNormalMag * sqrt_r2_2_minus_r2 * (num / den);
        // Real b_div_a = sqrt(den / (ny * ny + den));
        // More stable than computing b/a directly
        Real b_div_a = sqrt(den / (planeNormalMagSq * r2_2_minus_r2));

        // *Much* more numerically stable than computing a and c separately in
        // the (near) parabola cases. See MajorVertexStableVersion.nb
        Real a_plus_c  = planeIsolevel * (ny * sqrt_r2_2_minus_r2 + nx * (1 - c_r2)) - planeNormalMagSq;
             a_plus_c /= planeNormalMag * (nx * sqrt_r2_2_minus_r2 - ny * (1 - c_r2));

        // Real c = planeIsolevel * nx * ny / (den * planeNormalMag)
        //        - planeNormalMag * ny * oneMinusR2 / den;

        m_icurve = ConicSection<Real>(a, b_div_a, a_plus_c);
        m_a = a;
        m_c = a_plus_c - a;
    }

    // (t, z) coordinates of point on intersection plane closest to p
    std::pair<Real, Real> closestPlanePtCoords(const Pt &p) const {
        return { (p - m_iplaneOrigin).dot(m_dt),
                 (p - m_iplaneOrigin).dot(m_dz) };
    }

    Pt closestPlanePt(const Pt &p) const {
        return iplanePoint(closestPlanePtCoords(p));
    }

    Pt closestICurvePoint(Pt p) const {
        return iplanePoint(m_icurve.closestPtCoords(closestPlanePtCoords(p)));
    }

    // Evaluate the (t, z) parametrization of the intersection plane
    Pt iplanePoint(Real t, Real z) const {
        return m_iplaneOrigin + t * m_d_dt + z * m_d_dz;
    }
    Pt iplanePoint(const std::pair<Real, Real> &tz) const {
        return iplanePoint(tz.first, tz.second);
    }

    std::vector<Pt> intersectionCurvePoints(size_t nsubdiv, Real zSign = 1.0) const {
        std::vector<Pt> pts;
        // Ellipse spans from c - a to c + a on t axis.
        // generate points for t in max(-2, c - a), c + a
        Real tmin = std::max(-2.0, m_c - m_a), tmax = m_c + m_a;
        Real r1 = m_d_dt.norm();
        for (int i = 0; i < nsubdiv; ++i) {
            Real t = tmin + i * (tmax - tmin) / (nsubdiv - 1);
            Real z = m_icurve.z(t);
            if (isnan(z)) continue;
            pts.emplace_back(iplanePoint(t,  zSign * z));
        }
        return std::move(pts);
    }

    const Pt  &iplaneOrigin() const { return m_iplaneOrigin; }
    const Vec &iplaneNormal() const { return m_iplaneNormal; }
    // Intersection plane coordinate basis vectors.
    // These are the tangent vectors corresponding to points (1, 0) and (0, 1)
    // in the (t, z) coordinate system (relative to iplaneOrigin())
    const Vec &tBasisVec() const { return m_d_dt; }
    const Vec &zBasisVec() const { return m_d_dz; }

    const ConicSection<Real> &icurve() const { return m_icurve; }

private:
    Real m_c, m_a;
    bool m_degenerate = false;
    // Coordinate differentials for the intersection plane parametrized by
    // (t, z). Note that (t, z) coordinates are in a scaled space where the
    // joint radius is 1.
    Vec m_dt, m_dz;
    // Coordinate basis vectors tangent to the intersection plane.
    Vec m_d_dt, m_d_dz;

    // Point in 3D corresponding to intersection plane coordinates (0, 0)
    Pt m_iplaneOrigin;
    Vec m_iplaneNormal;

    // Intersection curve
    ConicSection<Real> m_icurve;

    // All intersection geometry formulas were derived in the special case where
    // edges have unit length and the joint radius is 1. This joint
    // configuration is then described by just three parameters: the angle
    // between edges and the two edges' radii at unit distance away. We call
    // these the edge's canonical parameters.
    // The canonical angle alpha is the same as the angle between edges in 3D.
    // However, we must compute the canonical endcap radii for both edges. See
    // EdgeIntersection.nb for more details.
    std::pair<Real, Real>
    m_canonicalParamsFromParams(Real l1, Real l2, Real r1, Real r2, Real r3) {
        return { 1.0 + (r2 - r1) / l1,
                 1.0 + (r3 - r1) / l2 };
    }
};

template<typename Real, bool visualize = false>
Real raycastSmoothingEval(const Point3D &p, const BinaryJoint<> &bjoint,
                          const SD::Primitives::InflatedEdge<Real> &edge1,
                          const SD::Primitives::InflatedEdge<Real> &edge2,
                          Real smoothness) {
    // Note: wasted work here--should combine to a single call that gets
    // both the point and the tangent
    auto e1p = edge1.closestFrustumPoint(p);
    auto e2p = edge2.closestFrustumPoint(p);

    auto e1t = edge1.frustumAxialTangent(p);
    auto e2t = edge2.frustumAxialTangent(p);

    const auto &origin = bjoint.iplaneOrigin();
    const auto &normal = bjoint.iplaneNormal();
    // Determine where the ray (e1p + alpha e1t) intersects the joint
    // edge intersection plane:
    // ((e1p + alpha e1t) - origin).normal = 0
    // (e1p - origin).normal = -alpha e1t.normal
    // alpha = (origin - e1p).normal / e1t.normal
    Real rayNormalComponent = e1t.dot(normal);
    // Frustum rays nearly parallel to the intersection plane should
    // never actually intersect the other edge
    Real modulation1, modulation2;
    Point3<Real> ipoint1, ipoint2;
    if (std::abs(rayNormalComponent) < 1e-5) { modulation1 = 0; }
    else {
        Real alpha = (origin - e1p).dot(normal) / rayNormalComponent;
        ipoint1 = e1p + e1t * alpha;
        modulation1 = tanh(2 * (1.0 - edge1.normal(ipoint1).dot(edge2.normal(ipoint1))));
    }

    rayNormalComponent = e2t.dot(normal);
    // Frustum rays nearly parallel to the intersection plane should
    // never actually intersect the other edge
    if (std::abs(rayNormalComponent) < 1e-5) { modulation2 = 0; }
    else {
        Real alpha = (origin - e2p).dot(normal) / rayNormalComponent;
        ipoint2 = e2p + e2t * alpha;
        modulation2 = tanh(2 * (1.0 - edge1.normal(ipoint2).dot(edge2.normal(ipoint2))));
    }

    if (std::abs(modulation1) < 1e-14) modulation1 = 0;
    if (std::abs(modulation2) < 1e-14) modulation2 = 0;
    if ((modulation1 < 0) || (modulation2 < 0)) {
        std::cout << "modulation1: " << modulation1 << std::endl;
        std::cout << "modulation2: " << modulation2 << std::endl;
        throw std::runtime_error("Negative modulation");
    }

    if (visualize) {
        drawArrow(e1p, e1p + 0.08 * e1t, 0.02, RGBColorf(1.0, 0.0, 0.0));
        drawArrow(e2p, e2p + 0.08 * e2t, 0.02, RGBColorf(0.0, 0.0, 1.0));

        drawArrow(ipoint1, ipoint1 + 0.08 * edge1.normal(ipoint1), 0.02, RGBColorf(0.0, 1.0, 0.0));
        drawArrow(ipoint1, ipoint1 + 0.08 * edge2.normal(ipoint1), 0.02, RGBColorf(0.0, 1.0, 0.0));

        drawArrow(ipoint2, ipoint2 + 0.08 * edge1.normal(ipoint2), 0.02, RGBColorf(1.0, 1.0, 0.0));
        drawArrow(ipoint2, ipoint2 + 0.08 * edge2.normal(ipoint2), 0.02, RGBColorf(1.0, 1.0, 0.0));
    }

    Real modulation = sqrt(modulation1 * modulation2);

    if (modulation * smoothness < 1/512.0) {
        return std::min(edge1.signedDistance(p), edge2.signedDistance(p));
    }
    return SD::exp_smin_reparam(edge1.signedDistance(p),
                                edge2.signedDistance(p),
                                modulation * smoothness);
}

// Interpolates from 1 at x = 0 to 0 at x = h
template<typename Real>
Real cubicSplineBump(Real x, Real h) {
    Real q = 2.0 * x / h;
    Real val = 0.0;

    if (q < 1.0)
        val = 1.0 + (.75 * q - 1.5) * q * q;
    else if (q < 2.0) {
        Real twomq = 2 - q;
        val = 0.25 * twomq * twomq * twomq;
    }
    return val;
}

// Interpolates from 1 at x = 0 to ~= 0 at x = h
template<typename Real>
Real expBump(Real x, Real h) {
    Real q = 2.0 * x / h;
    // return exp(-pow(std::abs(q), 2.5));
    return exp(-pow(std::abs(q), 4.0));
}

struct UnionedEdges {
    using Real = double;

    UnionedEdges(const GeometryParams &gp)
        : edge1(Point3D(0, 0, 0), Point3D(gp.edgeLen1, 0, 0), gp.radius1, gp.radius2), smoothingMode(gp.smoothingMode), smoothness(gp.smoothness) {
        hullModulationWidth = gp.hullModulationWidth;
        Point3D dir2(cos(gp.alpha), sin(gp.alpha), 0);
        dir2 *= gp.edgeLen2;
        edge2.set(Point3D(0, 0, 0), dir2, gp.radius1, gp.radius3);
        if ((smoothingMode == SMOOTH_RAYCAST) || (smoothingMode == SMOOTH_ELLIPSECP))
            bjoint.set(edge1.c1(), edge1.c2(), edge2.c2(), edge1.r1(), edge1.r2(), edge2.r2());
        convexHullInflatedTri.set(Point3D(0, 0, 0), Point3D(gp.edgeLen1, 0, 0), dir2,
                                  gp.radius1, gp.radius2, gp.radius3);
    }

    Real signedDistance(const Point3D &p) const {
        if (smoothingMode == SMOOTH_ELLIPSECP) {
            auto ccp = bjoint.closestICurvePoint(p);
            auto e1n = edge1.normal(ccp);
            auto e2n = edge2.normal(ccp);

            Real smoothingAmt = smoothness * tanh(2 * (1.0 - e1n.dot(e2n)));
            assert(!std::isnan(smoothingAmt));
            assert(!std::isinf(smoothingAmt));
            if (smoothingAmt < 1/512.0) {
                return std::min(edge1.signedDistance(p), edge2.signedDistance(p));
            }

            return SD::exp_smin_reparam(edge1.signedDistance(p),
                                        edge2.signedDistance(p),
                                        smoothingAmt);
        }
        else if (smoothingMode == SMOOTH_RAYCAST) {
            return raycastSmoothingEval(p, bjoint, edge1, edge2, smoothness);
        }
        else if (smoothingMode == SMOOTH_FULL) {
            return SD::exp_smin_reparam(edge1.signedDistance(p),
                                        edge2.signedDistance(p),
                                        smoothness);
        }
        else if (smoothingMode == SMOOTH_HULL) {
            return SD::exp_smin_reparam_accurate(edge1.signedDistance(p),
                                                 edge2.signedDistance(p),
                                                 smoothness * hullSmoothingModulation(p));
        }
        else if (smoothingMode == SMOOTH_NONE) {
            return std::min(edge1.signedDistance(p), edge2.signedDistance(p));
        }
        else throw std::runtime_error("Unknown smoothing mode");
    }

    Real hullSmoothingModulation(const Point3D &p) const {
        Real hullDist = convexHullInflatedTri.signedDistance(p);
        // Real modulation = 1.0 - cubicSplineBump(-std::min<Real>(hullDist, 0),
        //                                         hullModulationWidth);
        // Real modulation = 1.0 - expBump(std::min<Real>(hullDist, 0),
        //                                 hullModulationWidth);
        // if (hullDist < 0.0) {
            // Real z = std::max<Real>(1.0 + hullDist / edge1.r1(), 0.0); // from 0 at "center" to 1 at outside
            // z = sin(z * M_PI / 2.0);
            // return tanh(2 * (1.0 - z * z));
            Real z = 1.0 + (hullDist / edge1.r1()); // from 0 at "center" to 1 at outside
            z = std::max<Real>(z, 0.0);
            z *= 1.05;
            return 1.0 - tanh(pow(std::abs(z), 8.0));
        // }
        // No smoothing outside
        return 0.0;
    }

    void dumpModulationField(const std::string &path) const {
        std::vector<MeshIO::IOVertex > gv;
        std::vector<MeshIO::IOElement> ge;
        size_t gs = 50;
        std::vector<size_t> grid_sizes({gs, gs, gs});
        gen_grid(grid_sizes, gv, ge);
        ScalarField<Real> modulation(gv.size());
        for (size_t i = 0; i < gv.size(); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                gv[i].point[j] *= 2.0 / gs;
                gv[i].point[j] -= 1.0;
            }
            modulation[i] = hullSmoothingModulation(gv[i].point);
        }
        MSHFieldWriter writer(path, gv, ge);
        writer.addField("modulation", modulation, DomainType::PER_NODE);
    }

    SD::Primitives::InflatedEdge<Real> edge1, edge2;
    BBox<Point3D> boundingBox() const { return BBox<Point3D>(Point3D(-1, -1, -1), Point3D(1, 1, 1)); }
    SmoothingMode smoothingMode;
    Real smoothness, hullModulationWidth;
    BinaryJoint<> bjoint;
    SD::Primitives::InflatedTriangle<Real>  convexHullInflatedTri;
};

struct UnionedEdgeHull {
    using Real = double;

    UnionedEdgeHull(const GeometryParams &gp) {
        Point3D c1(0, 0, 0);
        Point3D c2(gp.edgeLen1, 0, 0);
        Point3D c3(cos(gp.alpha), sin(gp.alpha), 0);
        c3 *= gp.edgeLen2;
        inflatedTri.set(c1, c2, c3, gp.radius1, gp.radius2, gp.radius3);
    }

    Real signedDistance(const Point3D &p) const {
        return inflatedTri.signedDistance(p);
    }

    BBox<Point3D> boundingBox() const { return BBox<Point3D>(Point3D(-1, -1, -1), Point3D(1, 1, 1)); }
    SD::Primitives::InflatedTriangle<Real> inflatedTri;
};

// Shapes scale
float g_Zoom = 1.0f;
// Shape orientation (stored as a quaternion)
float g_Rotation[] = { 0.0f, 0.0f, 0.0f, 1.0f };
// Shapes material
float g_MatAmbient[] = { 0.4f, 0.4f, 0.4f, 1.4f };
float g_MatDiffuse[] = { 0.75f, 0.75f, 0.75f, 0.5f };
// Light parameter
float g_LightMultiplier = 1.0f;
float g_LightDirection[] = { 0.0f, 0.0f, -1.0f };
bool g_useVtxNormals = true;
bool g_liveRemesh = true;

bool g_showJoint = true;
bool g_showConvexHull = false;
bool g_showBlendingRegion = true;
bool g_showHullPlanes = false;
bool g_showIntersectionPlane = true;

float g_queryPtEllipseAngle = 0.30;
float g_queryPtDist = 0.23f;
float g_queryPtPerpDist = 0.0f;


struct TriMesh {
    TriMesh() { }
    TriMesh(const std::vector<MeshIO::IOVertex > &v,
            const std::vector<MeshIO::IOElement> &e) { setGeometry(v, e); }

    void setGeometry(const std::vector<MeshIO::IOVertex > &v,
                     const std::vector<MeshIO::IOElement> &e) {
        vertices = v;
        elements = e;
        recomputeNormals();
    }

    // Compute per-vertex, per-tri normals
    void recomputeNormals() {
        vtxNormals.assign(vertices.size(), Point3D::Zero());
        triNormals.clear(), triNormals.reserve(elements.size());
        for (const auto &e : elements) {
            Point3D v0 = vertices[e[0]].point,
                    v1 = vertices[e[1]].point,
                    v2 = vertices[e[2]].point;
            auto n = (v1 - v0).cross(v2 - v0).eval(); // area * normal
            for (size_t i = 0; i < 3; ++i)
                vtxNormals[e[i]] += n;
            triNormals.emplace_back(n / n.norm());
        }
        for (size_t i = 0; i < vtxNormals.size(); ++i)
            vtxNormals[i] /= vtxNormals[i].norm();
    }

    void clear() {
        vertices.clear(), elements.clear();
    }

    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;
    std::vector<Point3D> vtxNormals;
    std::vector<Point3D> triNormals;
};

TriMesh g_jointMesh, g_conicMesh1, g_conicMesh2, g_convexHullMesh, g_blendingRegionMesh;

// VCGSurfaceMesher<UnionedEdges> mesher;
IGLSurfaceMesherMC<UnionedEdges> mesher;
IGLSurfaceMesherMC<ConeQuadric<>> conicMesher;
IGLSurfaceMesherMC<UnionedEdgeHull> blendingRegionMesher;

GeometryParams g_geometryParams;
std::unique_ptr<GeometryParams> g_meshedGParams;

std::pair<ConeQuadric<>, ConeQuadric<>> getQuadrics() {
    return std::make_pair(
            ConeQuadric<>::FromXYEdge(0.0, g_geometryParams.edgeLen1, g_geometryParams.radius1, g_geometryParams.radius2),
            ConeQuadric<>::FromXYEdge(g_geometryParams.alpha, g_geometryParams.edgeLen2, g_geometryParams.radius1, g_geometryParams.radius3));
}

// Remesh only if the geometry has changed
void remesh() {
    if (g_meshedGParams && (*g_meshedGParams == g_geometryParams)) return;
    g_meshedGParams = Future::make_unique<GeometryParams>(g_geometryParams);

    const auto &gp = g_geometryParams;

    mesher.meshingOptions.marchingCubesGridSize = gp.mcGridSize;
    
    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;
    mesher.mesh(UnionedEdges(gp), vertices, elements);
    g_jointMesh.setGeometry(vertices, elements);

    // convexHull(g_jointMesh.vertices, g_convexHullMesh.vertices,
    //            g_convexHullMesh.elements);
    // g_convexHullMesh.recomputeNormals();

    blendingRegionMesher.meshingOptions = mesher.meshingOptions;
    blendingRegionMesher.mesh(UnionedEdgeHull(gp), vertices, elements, gp.isolevel);
    g_blendingRegionMesh.setGeometry(vertices, elements);

#if 0
    auto quadrics = getQuadrics();
    conicMesher.meshingOptions = mesher.meshingOptions;
    conicMesher.mesh(quadrics.first, vertices, elements);
    g_conicMesh1.setGeometry(vertices, elements);

    conicMesher.meshingOptions = mesher.meshingOptions;
    conicMesher.mesh(quadrics.second, vertices, elements);
    g_conicMesh2.setGeometry(vertices, elements);
#endif
}

// Remeshes only if the geometry has changed and live remeshing turned on
void liveRemesh() {
    if (g_meshedGParams && !g_liveRemesh) return;
    remesh();
}

const std::vector<size_t> &depthSortedOrder(
    const std::vector<MeshIO::IOVertex > &vertices,
    const std::vector<MeshIO::IOElement> &elements)
{
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat); // get column-major format

    // Sort by the z coordinate of the (unnormalized) homogenous vector
    // So we need only apply the third row of the matrix)
    Point3D projZ(mat[2], mat[6], mat[10]);
    double t = mat[14];

    static std::vector<Real> zcoord;
    zcoord.clear();
    zcoord.reserve(elements.size());
    for (const auto &e : elements) {
        Point3D barycenter = vertices[e[0]].point + vertices[e[1]].point + vertices[e[2]].point;
        // Real z = projZ.dot(barycenter) / 3.0 + t;
        Real z = projZ.dot(barycenter); // translation/uniform scaling doesn't change ordering
        zcoord.push_back(z);
    }

    static std::vector<size_t> perm;
    sortPermutation(zcoord, perm);
    return perm;
}


void drawArrow(const Point3D &tail, const Point3D &tip,
               Real radius, const RGBColorf &color) {
    Vector3D v = tip - tail;
    double len = v.norm();
    v /= len;
    // Get rotation mappping the z axis to v
    Eigen::Matrix<double, 4, 4> mat;
    mat.setIdentity();

    // z cross v
    Vector3D perp1(-v[1], v[0], 0.0);
    if (perp1.norm() < 1e-4) {
        // use y cross v instead if v ~= z (getting ~= x)
        perp1 = Vector3D(v[2], 0, -v[0]);
    }
    perp1 /= perp1.norm();
    Vector3D perp2(v.cross(perp1));
    mat.block<3, 1>(0, 0) = perp1;
    mat.block<3, 1>(0, 1) = perp2;
    mat.block<3, 1>(0, 2) = v;

    glPushMatrix();
    glTranslated(tail[0], tail[1], tail[2]);
    glMultMatrixd(mat.data());

    double arrowShaftRadius = 0.5 * radius;
    double arrowHeadLen = radius / tan(M_PI / 6.0);

    double shaftLen = len - arrowHeadLen;
    const size_t nsubdiv = 16;

    glColor3f(1, 0, 0);
    float matAmbient[] = {0.5f * color.r, 0.5f * color.g, 0.5f * color.b, 1.0f};
    float matDiffuse[] = {1.0f * color.r, 1.0f * color.g, 1.0f * color.b, 1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmbient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiffuse);

    // Draw shaft, if it exists
    if (shaftLen > 0) {
        GLUquadric *quad = gluNewQuadric();
        assert(quad);
        gluCylinder(quad,
                    arrowShaftRadius,
                    arrowShaftRadius,
                    shaftLen,
                    nsubdiv,
                    1);
        gluDeleteQuadric(quad);
    }
    glTranslated(0, 0, std::max<double>(shaftLen, 0));
    glutSolidCone(radius, arrowHeadLen, nsubdiv, 1);

    glPopMatrix();
}

void drawQueryPointInfo(const Point3D &queryPt) {
    glDisable(GL_LIGHTING);
    glColor3f(0.2, 1.0, 0.0);
    glPointSize(8.0);
    glBegin(GL_POINTS);
        glColor3f(0.2, 1.0, 0.0);
        glVertex3dv(queryPt.data());
    glEnd();
    glEnable(GL_LIGHTING);

    UnionedEdges joint(g_geometryParams);

    auto e1p = joint.edge1.closestPoint(queryPt);
    auto e2p = joint.edge2.closestPoint(queryPt);
    // auto e1p = joint.edge1.closestFrustumBorderPoint(queryPt, 0);
    // auto e2p = joint.edge2.closestFrustumBorderPoint(queryPt, 0);

    auto e1n = joint.edge1.normal(e1p);
    auto e2n = joint.edge2.normal(e2p);

    auto e1t = joint.edge1.frustumAxialTangent(e1p);
    auto e2t = joint.edge2.frustumAxialTangent(e2p);

    drawArrow(e1p, e1p + 0.25 * e1n, 0.04, RGBColorf(1.0, 0.0, 0.0));
    drawArrow(e2p, e2p + 0.25 * e2n, 0.04, RGBColorf(0.0, 0.0, 1.0));
}

void render(const TriMesh &m) {
    // Must come before glBegin: apparently cannot access modelview matrix after
    const std::vector<size_t> &order = depthSortedOrder(m.vertices, m.elements);

    glBegin(GL_TRIANGLES);
    for (size_t ei : order) {
        if (!g_useVtxNormals) glNormal3dv(m.triNormals[ei].data());
        for (size_t c : m.elements[ei]) {
            if (g_useVtxNormals) glNormal3dv(m.vtxNormals[c].data());
            glVertex3dv(m.vertices[c].point.data());
        }
    }
    glEnd();
}

void drawString(const std::string &str) {
    for (char c : str)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, c);
}

// Draw plane passing through o with tangent vectors t1 and t2.
void drawPlane(const Point3D &o, const Point3D &t1, const Point3D &t2) {
    glBegin(GL_QUADS);
        glVertex3dv((o - t1 - t2).eval().data());
        glVertex3dv((o + t1 - t2).eval().data());
        glVertex3dv((o + t1 + t2).eval().data());
        glVertex3dv((o - t1 + t2).eval().data());
    glEnd();
}

// Callback function called by GLUT to render screen
void Display(void) {
    float v[4]; // will be used to set light parameters
    float mat[4*4]; // rotation matrix

    // Clear frame buffer
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);

    // glShadeModel(GL_FLAT);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Set light
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    v[0] = v[1] = v[2] = g_LightMultiplier*0.4f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_AMBIENT, v);
    v[0] = v[1] = v[2] = g_LightMultiplier*0.8f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
    v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
    glLightfv(GL_LIGHT0, GL_POSITION, v);


    glPushMatrix();
    ConvertQuaternionToMatrix(g_Rotation, mat);
    glMultMatrixf(mat);
    glScalef(g_Zoom, g_Zoom, g_Zoom);


    // Opaque objects--must be drawn before translucent objects!
    if (g_showHullPlanes) {
        const auto &gp = g_geometryParams;
        Point3D c1(0.0, 0.0, 0.0),
                c2(gp.edgeLen1, 0, 0),
                c3(gp.edgeLen2 * Point3D(cos(gp.alpha), sin(gp.alpha), 0));

        // Determine the tangent plane for the three spheres
        Real nx = (gp.radius1 - gp.radius2) / gp.edgeLen1;
        Real ny = (gp.radius1 - gp.radius3 - c3[0] * nx) / c3[1];
        Real xyNormSq = nx * nx + ny * ny;
        if (xyNormSq <= 1.0) {
            Real nz = sqrt(1.0 - xyNormSq);

            // Draw the triangle formed by the sphere centers
            glDisable(GL_LIGHTING);
            glColor3f(0.0, 1.0, 0.0);
            glBegin(GL_TRIANGLES);
                glVertex3dv(c1.data());
                glVertex3dv(c2.data());
                glVertex3dv(c3.data());
            glEnd();

            Vector3D n(nx, ny, nz);
            Vector3D t1_candidate1 = n.cross(Vector3D(1, 0, 0));
            Vector3D t1_candidate2 = n.cross(Vector3D(0, 1, 0));
            Vector3D t1, t2;
            Real t1_candidate1_norm = t1_candidate1.norm();
            Real t1_candidate2_norm = t1_candidate2.norm();
            if (t1_candidate1_norm > t1_candidate2_norm)
                 t1 = t1_candidate1 / t1_candidate1_norm;
            else t1 = t1_candidate2 / t1_candidate2_norm;
            t2 = n.cross(t1);

            // Sphere tangent points
            // c1 + n * r1, ...
            Point3D pt1 = gp.radius1 * n;
            Point3D pt2 = c2 + gp.radius2 * n;
            Point3D pt3 = c3 + gp.radius3 * n;
            glColor3f(1.0, 1.0, 1.0);
            glPointSize(3.0);
            glBegin(GL_POINTS);
                glVertex3dv(pt1.data());
                glVertex3dv(pt2.data());
                glVertex3dv(pt3.data());
            glEnd();
            glColor3f(1.0, 0.0, 0.0);
            glBegin(GL_TRIANGLES);
                glVertex3dv(pt1.data());
                glVertex3dv(pt2.data());
                glVertex3dv(pt3.data());
            glEnd();

            glColor3f(1.0, 1.0, 1.0);
            glBegin(GL_LINES);
                glVertex3dv( c1.data());
                glVertex3dv(pt1.data());
                glVertex3dv( c2.data());
                glVertex3dv(pt2.data());
                glVertex3dv( c3.data());
                glVertex3dv(pt3.data());
            glEnd();


            // drawPlane((pt1 + pt2 + pt3) / 3.0, 0.25 * t1, 0.25 * t2);
            glEnable(GL_LIGHTING);
        }
    }

    // Draw the intersection plane
    if (g_showIntersectionPlane) {
        const auto &gp = g_geometryParams;
        BinaryJoint<> bjoint(Point3D(0, 0, 0),
                             Point3D(gp.edgeLen1, 0, 0),
                             gp.edgeLen2 * Point3D(cos(gp.alpha), sin(gp.alpha), 0),
                             gp.radius1, gp.radius2, gp.radius3);

        auto o = bjoint.iplaneOrigin();
        auto t = bjoint.tBasisVec();
        auto z = bjoint.zBasisVec();
        // t *= 0.5 / t.norm();
        // z *= 0.5 / z.norm();

        // {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}

        const auto &icurve = bjoint.icurve();
        auto ept_tz = icurve.pointAtAngle(g_queryPtEllipseAngle);
        auto enormal_tz = icurve.normalAtPoint(ept_tz);
        Point3D ellipsePt = bjoint.iplanePoint(ept_tz);
        Vector3D enormal = (t * enormal_tz.first + z * enormal_tz.second);
        enormal /= enormal.norm();
        Point3D queryPt = ellipsePt + g_queryPtDist * enormal;
        queryPt += g_queryPtPerpDist * z.cross(t) / (t.norm() * z.norm());

        auto cpp = bjoint.closestPlanePt(queryPt);
        auto ccp = bjoint.closestICurvePoint(queryPt);

        glDisable(GL_LIGHTING);
        glEnable(GL_POINT_SMOOTH);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        drawQueryPointInfo(queryPt);

        // glDisable(GL_LIGHTING);
        // glColor3f(1, 1, 1);
        // glPointSize(8.0);
        // glLineWidth(1.0);
        // glBegin(GL_POINTS);
        //     glVertex3dv(cpp.data());
        //     glVertex3dv(ccp.data());
        // glEnd();
        // glBegin(GL_LINES);
        //     glVertex3dv(ccp.data());
        //     glColor3f(0.2, 1.0, 0.0);
        //     glVertex3dv(cpp.data());

        //     glVertex3dv(cpp.data());
        //     glColor3f(0.2, 1.0, 0.0);
        //     glVertex3dv(queryPt.data());
        // glEnd();

        // // Draw the intersection curve (+/- z)
        // glColor3f(0.0, 0.2, 1.0);
        // glLineWidth(4);
        // glBegin(GL_LINE_STRIP);
        //     for (const auto &pt : bjoint.intersectionCurvePoints(1000, 1.0))
        //         glVertex3dv(pt.data());
        // glEnd();
        // glBegin(GL_LINE_STRIP);
        //     for (const auto &pt : bjoint.intersectionCurvePoints(1000, -1.0))
        //         glVertex3dv(pt.data());
        // glEnd();

        glColor3f(0.55, 0.1, 0.55);
        // drawPlane(o, 4 * t, 4 * z);
        glEnable(GL_LIGHTING);
        // drawArrow(o, o + 0.08 * bjoint.iplaneNormal(), 0.02, RGBColorf(1.0, 0.0, 0.0));


        UnionedEdges joint(g_geometryParams);
        auto e1n = joint.edge1.normal(queryPt);
        auto e2n = joint.edge2.normal(queryPt);

        // raycastSmoothingEval<double, true>(queryPt, bjoint, joint.edge1, joint.edge2, g_geometryParams.smoothness);

        glColor3f(1.0, 1.0, 1.0);
        glWindowPos2i(0, 0);
        drawString(std::to_string(joint.edge1.normal(queryPt).dot(joint.edge2.normal(queryPt))) + ", "
                 + std::to_string(joint.edge1.normal(ccp    ).dot(joint.edge2.normal(ccp    ))));

        // if (!g_meshedGParams || !(*g_meshedGParams == g_geometryParams)) {
        //     std::cout << o << std::endl;
        //     std::cout << t << std::endl;
        //     std::cout << z << std::endl;
        // }
    }

    // Translucent objects--must be depth sorted!
    // Geometry material
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, g_MatAmbient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, g_MatDiffuse);

    liveRemesh();

    if (g_showJoint)
        render(g_jointMesh);

    if (g_showBlendingRegion) {
        glEnable(GL_LIGHTING);
        RGBColorf color(0.5, 0.5, 0.5);
        float matAmbient1[] = {0.5f * color.r, 0.5f * color.g, 0.5f * color.b, 0.25f};
        float matDiffuse1[] = {1.0f * color.r, 1.0f * color.g, 1.0f * color.b, 0.25f};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmbient1);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiffuse1);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(-1.0, -1.0);
        render(g_blendingRegionMesh);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }

#if 0
    if (g_showConvexHull) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(-2.0, -2.0);

        glDisable(GL_LIGHTING);
        glLineWidth(1.0);
        glColor4f(1.0, 1.0, 1.0, 0.5);
        render(g_convexHullMesh);

        glDisable(GL_POLYGON_OFFSET_LINE);

        glEnable(GL_LIGHTING);
        RGBColorf color(0.2, 0.2, 1.0);
        float matAmbient1[] = {0.5f * color.r, 0.5f * color.g, 0.5f * color.b, 0.25f};
        float matDiffuse1[] = {1.0f * color.r, 1.0f * color.g, 1.0f * color.b, 0.25f};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmbient1);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiffuse1);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(-1.0, -1.0);
        render(g_convexHullMesh);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
#endif

#if 0
    RGBColorf color(0.2, 0.2, 1.0);
    float matAmbient1[] = {0.5f * color.r, 0.5f * color.g, 0.5f * color.b, 0.25f};
    float matDiffuse1[] = {1.0f * color.r, 1.0f * color.g, 1.0f * color.b, 0.25f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmbient1);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiffuse1);
    render(g_conicMesh1);

    color = RGBColorf(0.2, 1.0, 0.2);
    float matAmbient2[] = {0.5f * color.r, 0.5f * color.g, 0.5f * color.b, 0.25f};
    float matDiffuse2[] = {1.0f * color.r, 1.0f * color.g, 1.0f * color.b, 0.25f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmbient2);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiffuse2);
    render(g_conicMesh2);
#endif

    glPopMatrix();

    // Draw tweak bars
    TwDraw();

    // Present frame buffer
    glutSwapBuffers();
}


// Callback function called by GLUT when window size changes
void Reshape(int width, int height) {
    // Set OpenGL viewport and camera
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, (double)width/height, 1, 10);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.00, 0.00, 3.00,
              0.00, 0.00, 0.00,
              0.00, 1.00, 0.00);

    // Send the new window size to AntTweakBar
    TwWindowSize(width, height);

    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
    TwEventMouseButtonGLUT(button, state, x, y);
    glutPostRedisplay();
}

void MotionFunc(int x, int y) {
    bool handled = TwEventMouseMotionGLUT(x, y);
    if (handled) glutPostRedisplay();
}

void SpecialKeyboardFunc(int k, int x, int y) {
    TwEventSpecialGLUT(k, x, y);
    glutPostRedisplay();
}

void KeyboardFunc(unsigned char k, int x, int y) {
    TwEventKeyboardGLUT(k, x, y);
    glutPostRedisplay();
}

// Function called at exit
void Terminate(void) {
    TwTerminate();
}

void TW_CALL RemeshCallback(void *caller) {
    remesh();
}

void TW_CALL SaveCallback(void *caller) {
    MeshIO::save("geometry.msh", g_jointMesh.vertices, g_jointMesh.elements);
}

void TW_CALL SaveHullCallback(void *caller) {
    MeshIO::save("hull.msh", g_convexHullMesh.vertices, g_convexHullMesh.elements);
}

void TW_CALL SaveModulationCallback(void *caller) {
    UnionedEdges ue(g_geometryParams);
    ue.dumpModulationField("modulation.msh");
}

void TW_CALL SaveEllipseSD(void *caller) {
    const auto &gp = g_geometryParams;
    BinaryJoint<> bjoint(Point3D(0, 0, 0),
                         Point3D(gp.edgeLen1, 0, 0),
                         gp.edgeLen2 * Point3D(cos(gp.alpha), sin(gp.alpha), 0),
                         gp.radius1, gp.radius2, gp.radius3);
    const auto curvePts = bjoint.intersectionCurvePoints(1000, 1.0);
    UnionedEdges edges(g_geometryParams);
    auto quadrics = getQuadrics();

    std::ofstream ellipseSD("ellipse_sd.txt");
    ellipseSD << std::setprecision(19);
    for (size_t i = 0; i < curvePts.size(); ++i) {
        const auto &p = curvePts[i];
        ellipseSD << i
                  << '\t' << edges.signedDistance(p)
                  << '\t' << edges.edge1.signedDistance(p)
                  << '\t' << edges.edge2.signedDistance(p)
                  << '\t' << quadrics.first.signedDistance(p)
                  << '\t' << quadrics.second.signedDistance(p)
                  << std::endl;
    }
}

int main(int argc, char *argv[]) {
    TwBar *bar; // Pointer to the tweak bar

    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1024, 768);
    glutCreateWindow("Inflated edge joint viewer");
    glutCreateMenu(NULL);

    // Set GLUT callbacks
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    atexit(Terminate);  // Called after glutMainLoop ends

    // Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);

    // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
    glutMouseFunc(MouseFunc);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc(MotionFunc);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
    glutPassiveMotionFunc(MotionFunc);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc(KeyboardFunc);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc(SpecialKeyboardFunc);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);

    // Create a tweak bar
    bar = TwNewBar("Settings");
    TwDefine("Settings size='200 400' color='96 96 96' "); // change default tweak bar size and color

    // Add 'g_Zoom' to 'bar': this is a modifable (RW) variable of type TW_TYPE_FLOAT. Its key shortcuts are [z] and [Z].
    TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &g_Zoom,
               " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

    // Add 'g_Rotation' to 'bar': this is a variable of type TW_TYPE_QUAT4F which defines the object's orientation
    TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &g_Rotation,
               " label='Object rotation' opened=true help='Change the object orientation.' ");

    // Add 'g_LightMultiplier' to 'bar': this is a variable of type TW_TYPE_FLOAT. Its key shortcuts are [+] and [-].
    TwAddVarRW(bar, "Multiplier", TW_TYPE_FLOAT, &g_LightMultiplier,
               " label='Light booster' min=0.1 max=4 step=0.02 keyIncr='+' keyDecr='-' help='Increase/decrease the light power.' group='Lighting'");
    TwAddVarRW(bar, "LightDir", TW_TYPE_DIR3F, &g_LightDirection,
               " label='Light direction' help='Change the light direction.' group='Lighting'");
    TwAddVarRW(bar, "Smooth shading", TW_TYPE_BOOLCPP, &g_useVtxNormals, "label='Smooth shading' group='Geometry' group='Lighting'");
    TwDefine("Settings/Lighting opened=false");
    TwAddVarRW(bar, "Ambient", TW_TYPE_COLOR3F, &g_MatAmbient, " group='Material' ");
    TwAddVarRW(bar, "Diffuse", TW_TYPE_COLOR3F, &g_MatDiffuse, " group='Material' ");

    TwAddVarRW(bar, "Ellipse Angle", TW_TYPE_FLOAT, &g_queryPtEllipseAngle, "label='Angle of ellipse closest point' min=0.0 max=6.28 step=0.02 group='Query point'");
    TwAddVarRW(bar, "Distance",  TW_TYPE_FLOAT, &g_queryPtDist, "label='Distance' min=-1 max=1 step=0.01 group='Query point'");
    TwAddVarRW(bar, "Perp Distance",  TW_TYPE_FLOAT, &g_queryPtPerpDist, "label='Distance perpendicular to intersection plane' min=-1 max=1 step=0.01 group='Query point'");

    TwAddVarRW(bar, "alpha", TW_TYPE_FLOAT, &(g_geometryParams.alpha),
               " label='2nd Edge Orientation' min=0.00 max=3.1415926 step=0.02 help='Change orientation angle of second edge' group='Geometry'");
    TwAddVarRW(bar, "radius1", TW_TYPE_FLOAT, &(g_geometryParams.radius1), " label='center radius'     min=0.05 max=0.5 step=0.01 group='Geometry' ");
    TwAddVarRW(bar, "radius2", TW_TYPE_FLOAT, &(g_geometryParams.radius2), " label='edge 1 tip radius' min=0.05 max=0.5 step=0.01 group='Geometry' ");
    TwAddVarRW(bar, "radius3", TW_TYPE_FLOAT, &(g_geometryParams.radius3), " label='edge 2 tip radius' min=0.05 max=0.5 step=0.01 group='Geometry' ");
    TwAddVarRW(bar, "edgeLen1", TW_TYPE_FLOAT, &(g_geometryParams.edgeLen1), " label='edge 1 length' min=0.05 max=1.0 step=0.01 group='Geometry' ");
    TwAddVarRW(bar, "edgeLen2", TW_TYPE_FLOAT, &(g_geometryParams.edgeLen2), " label='edge 2 length' min=0.05 max=1.0 step=0.01 group='Geometry' ");

    TwEnumVal smoohingModeEV[] = { {SMOOTH_NONE, "None"}, {SMOOTH_FULL, "Full"}, {SMOOTH_ELLIPSECP, "Ellipse CP"}, {SMOOTH_RAYCAST, "Raycast"}, {SMOOTH_HULL, "Hull"} };
    TwType smoothingModeTwEnum = TwDefineEnum("SmoothingMode", smoohingModeEV, 5);
    // Adding season to bar
    TwAddVarRW(bar, "Smoothing mode", smoothingModeTwEnum, &g_geometryParams.smoothingMode, "label='Joint smoothing mode' group='Geometry' ");
    TwAddVarRW(bar, "smoothness", TW_TYPE_FLOAT, &(g_geometryParams.smoothness), " label='Smoothing amount' min=0.001 max=0.2 step=0.01 group='Geometry' ");
    TwAddVarRW(bar, "Hull modulation width", TW_TYPE_FLOAT, &(g_geometryParams.hullModulationWidth), " label='Hull modulation width' min=0.001 max=0.5 step=0.001 group='Geometry' ");

    TwAddVarRW(bar, "MC Grid Size", TW_TYPE_INT32, &(g_geometryParams.mcGridSize), " min=32 max=1024 help='The marching cubes resolution' group='Meshing' ");
    TwAddVarRW(bar, "Live Remesh", TW_TYPE_BOOLCPP, &g_liveRemesh,
                   " group='Meshing' ");
    TwAddVarRW(bar, "Show joint mesh", TW_TYPE_BOOLCPP, &g_showJoint, " group='Meshing' ");
    TwAddVarRW(bar, "Show intersection plane", TW_TYPE_BOOLCPP, &g_showIntersectionPlane,
                   " group='Meshing' ");
    TwAddVarRW(bar, "Show blending region", TW_TYPE_BOOLCPP, &g_showBlendingRegion, " group='Meshing' ");
    TwAddVarRW(bar, "Show conve hull planes", TW_TYPE_BOOLCPP, &g_showHullPlanes, " group='Meshing' ");
    TwAddVarRW(bar, "Blending region isolevel", TW_TYPE_FLOAT, &g_geometryParams.isolevel, "label='Blending region isolevel' min=-0.75 max=0.75 step=0.01 group='Meshing'");
    TwAddVarRW(bar, "Show convex hull", TW_TYPE_BOOLCPP, &g_showConvexHull, " group='Meshing' ");
    TwAddButton(bar, "Remesh", RemeshCallback, NULL, "group='Meshing'");
    TwAddButton(bar, "Save geometry", SaveCallback, NULL, "group='Output'");
    TwAddButton(bar, "Save convex hull", SaveHullCallback, NULL, "group='Output'");
    TwAddButton(bar, "Save modulation field", SaveModulationCallback, NULL, "group='Output'");
    TwAddButton(bar, "Write ellipse signed distances", SaveEllipseSD, NULL, "group='Output'");

    // Call the GLUT main loop
    glutMainLoop();

    return 0;
}

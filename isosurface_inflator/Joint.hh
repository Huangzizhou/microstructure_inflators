////////////////////////////////////////////////////////////////////////////////
// Joint.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Represents and computes signed distance to the smoothed joint geometry
//      at a particular vertex.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/08/2016 23:55:55
////////////////////////////////////////////////////////////////////////////////
#ifndef JOINT_HH
#define JOINT_HH

#include "SphereConvexHull.hh"
#include <Future.hh>
#include "AutomaticDifferentiation.hh"

#include <stdexcept>

enum class JointBlendMode { FULL, HULL, HULL_HALF_EDGE };
template<typename Real>
class Joint {
public:
    Joint(const std::vector<Point3<Real>> &centers,
          const std::vector<Real>         &radii,
          Real blendingAmt, JointBlendMode mode)
    { setParameters(centers, radii, blendingAmt, mode); }

    // Assumes the first (center, radius) specifies the joint sphere
    void setParameters(std::vector<Point3<Real>> centers, // modified inside
                       std::vector<Real>         radii,   // modified inside
                       Real blendingAmt,
                       JointBlendMode blendMode) {
        m_mode = blendMode;
        assert(centers.size() == radii.size());
        if (centers.size() < 3) {
            {
                std::cout << "Joint constructor called on centers, radii:" << std::endl;
                for (const auto &pt : centers) {
                    std::cout << "{"
                        << pt[0] << ", "
                        << pt[1] << ", "
                        << pt[2] << "}, ";
                }
                std::cout << std::endl;
                std::cout << "{";
                for (const Real &r : radii)
                    std::cout << r << ", ";
                std::cout << "}" << std::endl;
            }
            throw std::runtime_error("Joint must comprise at least 3 spheres");
        }
        m_r1 = radii[0];
        m_c1 = centers[0];
        if (m_mode == JointBlendMode::HULL_HALF_EDGE) {
            // In HULL_HALF_EDGE mode, shrink the hull region to extend only
            // halfway through each edge so that neighboring joint blending
            // regions do not overlap. This should prevent creases from forming
            // when joint geometries are unioned together.
            for (size_t i = 1; i < centers.size(); ++i) {
                centers[i] += centers[0];
                radii[i]   += radii[0];
                centers[i] *= 0.5;
                radii[i]   *= 0.5;
            }
        }
        m_blendingAmt = blendingAmt;
        m_blendingHull = Future::make_unique<SD::Primitives::SphereConvexHull<Real>>(centers, radii);
    }

    template<typename Real2>
    Real2 smoothingAmt(const Point3<Real2> &p) const {
        if (m_mode == JointBlendMode::FULL) { return m_blendingAmt; }

        Real2 hullDist = m_blendingHull->signedDistance(p);
        Real2 z = 1.0 + (hullDist / m_r1); // from 0 at "center" to 1 at outside
        z = 1.05 * std::max<Real2>(z, 0.0);
        Real2 modulation = 1.0 - tanh(pow(z, 8.0));

        return modulation * m_blendingAmt;
    }

    Real blendParam() const { return m_blendingAmt; }

    Real         r1() const { return m_r1; }
    Point3<Real> c1() const { return m_c1; }

    SD::Primitives::SphereConvexHull<Real> blendingHull() const { assert(m_blendingHull); return *m_blendingHull; }

private:
    std::unique_ptr<SD::Primitives::SphereConvexHull<Real>> m_blendingHull;
    Real m_blendingAmt;
    // Center and radius of the joint sphere.
    Point3<Real> m_c1;
    Real         m_r1;
    JointBlendMode m_mode; 
};

#endif /* end of include guard: JOINT_HH */

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

#include <stdexcept>

template<typename Real>
class Joint {
public:
    enum class BlendMode { FULL, HULL };
    Joint(const std::vector<Point3<Real>> &centers,
          const std::vector<Real>         &radii,
          Real blendingAmt)
    { setParameters(centers, radii, blendingAmt); }

    // Assumes the first (center, radius) specifies the joint sphere
    void setParameters(const std::vector<Point3<Real>> &centers,
                       const std::vector<Real>         &radii,
                       Real blendingAmt) {
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
        m_blendingAmt = blendingAmt;
        m_blendingHull = Future::make_unique<SD::Primitives::SphereConvexHull<Real>>(centers, radii);
    }

    template<typename Real2>
    Real2 smoothingAmt(const Point3<Real2> &p) const {
        Real2 hullDist = m_blendingHull->signedDistance(p);
        Real2 z = 1.0 + (hullDist / m_r1); // from 0 at "center" to 1 at outside
        z = 1.05 * std::max<Real2>(z, 0.0);
        Real2 modulation = 1.0 - tanh(pow(z, 8.0));

        if (m_mode == BlendMode::FULL) modulation = 1.0;
        return modulation * m_blendingAmt;
    }

private:
    std::unique_ptr<SD::Primitives::SphereConvexHull<Real>> m_blendingHull;
    Real m_blendingAmt;
    Real m_r1;
    BlendMode m_mode = BlendMode::HULL; 
};

#endif /* end of include guard: JOINT_HH */

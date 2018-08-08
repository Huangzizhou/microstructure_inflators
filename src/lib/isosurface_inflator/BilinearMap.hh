////////////////////////////////////////////////////////////////////////////////
// BilinearMap.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Simple class for representing a bilinear map from a reference square [-1,1]²
//  to an arbitrary quadrilateral.
*/
//  Author:  Jérémie Dumas (jdumas), jeremie.dumas@ens-lyon.org
//  Created:  08/08/2018 11:00:00
////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "InflatorTypes.hh"

class BilinearMap {

public:
    BilinearMap() = default;

    template<typename T>
    BilinearMap(T pts) {
        a[0] = pts[0][0]; a[1] = pts[0][1];
        b[0] = pts[1][0]; b[1] = pts[1][1];
        c[0] = pts[2][0]; c[1] = pts[2][1];
        d[0] = pts[3][0]; d[1] = pts[3][1];
        a[2] = b[2] = c[2] = d[2] = 0.0;
    }

    template<typename Real>
    Point3<Real> apply(Real u, Real v) const {
        return a*(4*u*v - 4*u - 4*v + 4) + b*(-4*u*v + 4*u + 2*v - 2) + c*(4*u*v - 2*u - 2*v + 1) + d*(-4*u*v + 2*u + 4*v - 2);
    }

    Eigen::Matrix3d jacobian(double u, double v) const {
        Point3d dfdu = -2*a + 2*b + 2*(2*v - 1)*(a - b + c - d);
        Point3d dfdv = -2*a + 2*d + 2*(2*u - 1)*(a - b + c - d);
        Eigen::Matrix3d jac;
        jac <<
            dfdu[0], dfdv[0], 0,
            dfdu[1], dfdv[1], 0,
            dfdu[2], dfdv[2], 1;
        return jac;
    }

private:
    Point3d a, b, c, d;
};

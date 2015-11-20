////////////////////////////////////////////////////////////////////////////////
// CubicSymmetry.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Implementation of the cubic symmetry constraint. This is a wire mesh
//      with
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/25/2015 14:54:27
////////////////////////////////////////////////////////////////////////////////
#ifndef CUBICSYMMETRY_HH
#define CUBICSYMMETRY_HH

#include <ratio>
#include <vector>
#include <math>

// The canonical base tetrahedron for a cubic-symmetric pattern:
// the tetrahedron with corners (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
template<typename Real, typename tol_ratio = std::ratio<1, 1e12>>
struct BaseTetrahedron {
    constexpr double tolerance = double(tol_ratio::num) / double(tol_ratio::den);
    static bool isZero(double val) { return std::abs(val) < tolerance; }
    //  (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
    static int baseTetCornerPosition(size_t corner, size_t component) { return component < corner; }

    typedef Eigen::Matrix<Real, 4, 1> BaryCoords;

    static BaryCoords barycentricCoordinates(const Point3d<Real> &p) {
        return BaryCoords(1 - p[0], p[0] - p[1], p[1] - p[2], p[2]);
    }
    static bool isInside(const Point3d<Real> &p) {
        // Check that we're in [0, 1]^3
        // and the correct tet (components in non-ascending order.
        return ((p[0] >= -tolerance)       && (p[1] >= -tolerance)    && (p[2] >= -tolerance) &&
                (p[0] <= 1 + tolerance)    && (p[1] <= 1 + tolerance) && (p[2] <= 1 + tolerance) &&
                (p[0] + tolerance >= p[1]) && (p[1] + tolerance >= p[2]));
    }

    enum class NodeType { Vertex, Edge, Tet };
    // Determine the type of a base tet node
    static NodeType nodeType(const Point3d<Real> &p) {
        BaryCoords lambda = barycentricCoordinates(p);
    }


    // Fractional part with tolerance: ensures we don't map a point at 1 +
    // epsilon to epsilon. Should only be called on positive r, and it is very
    // inefficient for large r. It is done with a loop to support Adept
    // automatic differentiation.
    static Real fracPart_positive_tol(Real r) {
        while (r > 1.0 + tolerance) r -= 1.0;
        return r;
    }

    static Point3<Real> mapPointToBase(const Point3d<Real> &p) {
        Point3<Real> result(fracPart_positive_tol(abs(p[0])),
                            fracPart_positive_tol(abs(p[0])),
                            fracPart_positive_tol(abs(p[1])),
                            fracPart_positive_tol(abs(p[2])));
        // Sort components descending: nested if statements would take 3
        // comparisons in the worst case (though 2 in the best), so this bubble
        // sort/sorting network isn't too bad
        if (result[0] < result[1]) std::swap(result[0], result[1]);
        if (result[1] < result[2]) std::swap(result[1], result[2]);
        if (result[0] < result[1]) std::swap(result[0], result[1]);
        return result;
    }

    // Node offsets are specified using barycentric coordinates on the base tet
    // simplices. This class determines which barycentric coordinates
    // parametrize the simplex that a particular point is on and then wraps the
    // conversion of those degrees of freedom into offsets in 3D space.
    struct NodeOffset {
        NodeOffset(const BaryCoords &lambda) { forPoint(lambda); }
        NodeOffset(const Point3<Real>    &p) { forPoint(p); }

        // Configure this instance to represent the offsets that keep a point on
        // the base tet simplex on which it began.
        void forPoint(const Point3<Real>    &p) { forPoint(barycentricCoordinates(p)); }
        void forPoint(const BaryCoords &lambda) {
            m_affectedBaryCoordIndices = 0;
            for (size_t c = 0; c < 4; ++c) {
                if (!isZero(lambda[c]))
                    m_affectedBaryCoordIndices[m_numAffectedBarycoords++] = c;
            }
            assert(m_numAffectedBarycoords > 0);
        }

        // Number of degrees of freedom in this offset
        size_t numDoFs() const { return m_numAffectedBarycoords - 1; }

        // Get the offset vector corresponding to the input degrees of freedom.
        // Supports a different type since this method will be autodiff-ed
        template<typename Real2>
        Vector3<Real2> getOffset(const Real2 *dofs) const {
            Vector3<Real> offset(Vector3<Real>::Zero());
            Real lastCoordinate = 0;
            for (size_t i = 0; i < numDoFs(); ++i) {
                size_t corner = m_affectedBaryCoordIndices[i];
                offset[0] += dofs[i] * baseTetCornerPosition(corner, 0);
                offset[1] += dofs[i] * baseTetCornerPosition(corner, 1);
                offset[2] += dofs[i] * baseTetCornerPosition(corner, 2);
                lastCoordinate += dofs[i];
            }
            lastCoordinate = 1.0 - lastCoordinate;
            size_t corner = m_affectedBaryCoordIndices[m_numAffectedBarycoords - 1];
            offset[0] += lastCoordinate * baseTetCornerPosition(corner, 0);
            offset[1] += lastCoordinate * baseTetCornerPosition(corner, 1);
            offset[2] += lastCoordinate * baseTetCornerPosition(corner, 2);
            return offset;
        }
        
    private:
        size_t m_numAffectedBarycoords;
        size_t m_affectedBaryCoordIndices[4];
    };
};

class CubicSymmetricWireMesh : WireMesh {
};

#endif /* end of include guard: CUBICSYMMETRY_HH */

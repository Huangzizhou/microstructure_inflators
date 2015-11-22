////////////////////////////////////////////////////////////////////////////////
// NodePositioners.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Classes wrapping nodes' positional degrees of freedom in the symmetry
//      base unit. DoFs are constrained to keep nodes inside the base unit.
//
//      Positioners for two base units types are implemented: box (for triply
//      periodic and orthotropic symmetry), and tet (for cubic symmetry).
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/20/2015 11:44:59
////////////////////////////////////////////////////////////////////////////////
#ifndef NODEPOSITIONERS_HH
#define NODEPOSITIONERS_HH

#include <Geometry.hh>
#include "InflatorTypes.hh"
#include "FuzzySign.hh"

// Offsets must keep a point on the base cube region (face, edge, node, cube
// interior) on which it started.
// Position variables are always in the range [0, 1]
template<typename Real, typename TOL>
struct BoxNodePositioner {
    BoxNodePositioner(const BBox<Point3<Real>> &cell, const Point3<Real> &p) { forPoint(cell, p); }

    void forPoint(const BBox<Point3<Real>> &cell, const Point3<Real> &p) {
        for (size_t c = 0; c < 3; ++c) {
            // Node can move in any coordinate it doesn't share with a
            // min/max cell face.
            if      (isZero<TOL>(p[c] - cell.minCorner[c])) m_ctype[c] = ComponentType::MinFace;
            else if (isZero<TOL>(p[c] - cell.maxCorner[c])) m_ctype[c] = ComponentType::MaxFace;
            else                                            m_ctype[c] = ComponentType::Free;
        }
        m_dimensions = cell.dimensions();
        m_minCorner = cell.minCorner;
    }

    size_t numDoFs() const { return (m_ctype[0] == ComponentType::Free) +
                                    (m_ctype[1] == ComponentType::Free) +
                                    (m_ctype[2] == ComponentType::Free); }

    // Get the 3D position corresponding to the input degrees of freedom.
    // Supports a different type since this method will be autodiff-ed
    template<typename Real2>
    Point3<Real2> getPosition(const Real2 *dofs) const {
        Vector3<Real2> pos;
        size_t d = 0;
        for (size_t i = 0; i < 3; ++i) {
            pos[i] = m_minCorner[i];
            if (m_ctype[i] == ComponentType::Free)    pos[i] += dofs[d++] * m_dimensions[i];
            if (m_ctype[i] == ComponentType::MaxFace) pos[i] += m_dimensions[i];
        }
        assert(d == numDoFs());
        return pos;
    }

    // Get the degrees of freedom corresponding to a 3D position. Assumes that
    // the position lies within the space spanned by this positioner.
    template<typename Real2>
    void getDoFsForPoint(const Point3<Real> &p, Real2 *dofs) const {
        size_t d = 0;
        for (size_t i = 0; i < 3; ++i) {
            if (m_ctype[i] == ComponentType::Free)
                dofs[d++] = (p[i] - m_minCorner[i]) / m_dimensions[i];
            else {
                // Verify that p lies in the expected subspace
                if (m_ctype[i] == ComponentType::MinFace)
                    assert(isZero<TOL>(p[i] - m_minCorner[i]));
                if (m_ctype[i] == ComponentType::MaxFace)
                    assert(isZero<TOL>(p[i] - m_minCorner[i] - m_dimensions[i]));
            }
            
        }
        assert(d == numDoFs());
    }

private:
    enum class ComponentType { MinFace, MaxFace, Free };
    ComponentType m_ctype[3];
    Point3<Real> m_minCorner;
    Vector3<Real> m_dimensions;
};

// Offsets must keep a point on the base tetrahedron simplex on which it
// started (face, node, edge, tet interior). The positional degrees of
// freedom are barycentric coordinates for this simplex.
// Position variables are always in the range [0, 1], but there is an
// additional linear inequality constraint restricting their sum to 1.
template<typename Real, typename TOL>
struct TetNodePositioner {
    typedef Eigen::Matrix<Real, 4, 1> BaryCoords;
    //  (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
    static int baseTetCornerPosition(size_t corner, size_t component) { return component < corner; }
    static BaryCoords barycentricCoordinates(const Point3<Real> &p) {
        return BaryCoords(1 - p[0], p[0] - p[1], p[1] - p[2], p[2]);
    }

    TetNodePositioner(const BBox<Point3<Real>> & /* cell */, const Point3<Real> &p) { forPoint(p); }
    TetNodePositioner(const Point3<Real> &p) { forPoint(p); }

    // Assumes point is within the base unit!
    void forPoint(const Point3<Real>    &p) { forPoint(barycentricCoordinates(p)); }
    void forPoint(const BaryCoords &lambda) {
        m_numAffectedBarycoords = 0;
        for (size_t c = 0; c < 4; ++c) {
            if (!isZero<TOL>(lambda[c]))
                m_affectedBaryCoordIndices[m_numAffectedBarycoords++] = c;
        }
        assert(m_numAffectedBarycoords > 0);
    }

    size_t numDoFs() const { return m_numAffectedBarycoords - 1; }

    // Get the 3D position corresponding to the input degrees of freedom.
    // Supports a different type since this method will be autodiff-ed
    template<typename Real2>
    Point3<Real2> getPosition(const Real2 *dofs) const {
        Vector3<Real2> pos(Vector3<Real2>::Zero());
        Real2 lastCoordinate = 0;
        for (size_t i = 0; i < numDoFs(); ++i) {
            size_t corner = m_affectedBaryCoordIndices[i];
            pos[0] += dofs[i] * baseTetCornerPosition(corner, 0);
            pos[1] += dofs[i] * baseTetCornerPosition(corner, 1);
            pos[2] += dofs[i] * baseTetCornerPosition(corner, 2);
            lastCoordinate += dofs[i];
        }
        lastCoordinate = 1.0 - lastCoordinate;
        size_t corner = m_affectedBaryCoordIndices[m_numAffectedBarycoords - 1];
        pos[0] += lastCoordinate * baseTetCornerPosition(corner, 0);
        pos[1] += lastCoordinate * baseTetCornerPosition(corner, 1);
        pos[2] += lastCoordinate * baseTetCornerPosition(corner, 2);
        return pos;
    }

    // Get the degrees of freedom corresponding to a 3D position. Assumes that
    // the position lies within the space spanned by this positioner.
    template<typename Real2>
    void getDoFsForPoint(const Point3<Real> &p, Real2 *dofs) const {
        auto b = barycentricCoordinates(p);
        // Verify that b is in the base tet.
        for (size_t i = 0; i < 4; ++i)
            assert(isPositive<TOL>(b[i])), assert(isNegative<TOL>(b[i] - 1));
        Real simplexWeight = 0.0;
        for (size_t i = 0; i < m_numAffectedBarycoords; ++i) {
            simplexWeight += b[m_affectedBaryCoordIndices[i]];
            if (i < numDoFs()) dofs[i] = b[m_affectedBaryCoordIndices[i]];
        }
        // Verify that the point is within the expected simplex
        assert(isZero<TOL>(simplexWeight - 1.0));
    }

private:
    size_t m_numAffectedBarycoords;
    size_t m_affectedBaryCoordIndices[4];
};

#endif /* end of include guard: NODEPOSITIONERS_HH */

////////////////////////////////////////////////////////////////////////////////
// Symmetry.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Holds symmetry classes to be used as template parameters for the
//      WireMesh class. These provide the queries needed to determine the base
//      unit for a particular tiling symmetry and the degrees of freedom
//      parametrizing the symmetry-preserving subspace of patterns. They can
//      also generate (a finite subset of) the elements of the corresponding
//      symmetry group.
//
//      All symmetry queries are performed with a tolerance that is configurable
//      by template parameter.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/26/2015 17:55:53
////////////////////////////////////////////////////////////////////////////////
#ifndef SYMMETRY_HH
#define SYMMETRY_HH
#include <ratio>
#include <vector>
#include <memory>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <Geometry.hh>
#include "InflatorTypes.hh"
#include "Isometries.hh"
#include "NodePositioners.hh"
#include "FuzzySign.hh"
#include "AutomaticDifferentiation.hh"

namespace Symmetry {

// Enum value is the node's number of positional degrees of freedom.
enum class NodeType : unsigned int { Vertex = 0, Edge = 1, Face = 2, Interior = 3 };

// Forward declarations of Symmetry types.
template<typename TOL = DEFAULT_TOL> struct TriplyPeriodic;
template<typename TOL = DEFAULT_TOL> struct Orthotropic;
template<typename TOL = DEFAULT_TOL> struct Cubic;

// We need a traits class for CRTP to look up the correct NodePositioner class.
// This traits class must be specialized for each symmetry type.
// TriplyPeriodic and Orthotropic symmetries have a box base cell,
// Cubic symmetry has a tet base cell.
template<class Sym> struct SymmetryTraits { };
template<typename TOL> struct SymmetryTraits<TriplyPeriodic<TOL>> { template<typename Real> using NodePositioner = BoxNodePositioner<Real, TOL>; };
template<typename TOL> struct SymmetryTraits<Orthotropic<TOL>>    { template<typename Real> using NodePositioner = BoxNodePositioner<Real, TOL>; };
template<typename TOL> struct SymmetryTraits<Cubic<TOL>>          { template<typename Real> using NodePositioner = TetNodePositioner<Real, TOL>; };

// Implements some of the shared interface of the symmetry classes
template<class Sym>
struct SymmetryCRTP {
    template<typename Real>
    using NodePositioner = typename SymmetryTraits<Sym>::template NodePositioner<Real>;

    template<typename Real>
    static NodePositioner<Real> nodePositioner(const Point3<Real> &p) {
        assert(Sym::inBaseUnit(p));
        return NodePositioner<Real>(Sym::template representativeMeshCell<Real>(), p);
    }

    template<typename Real>
    static NodeType nodeType(const Point3<Real> &p) {
        return static_cast<NodeType>(nodePositioner(p).numDoFs());
    }
};

template<typename T, bool _isAutoDiffType = IsAutoDiffType<T>::value>
struct OptionalFMod2 {
    static void run(T &val) {
        T q = (int) (val / 2); // round toward zero
        // if (abs(q) > 1000) {
        //     std::cout << "queried far point " << p.transpose() << std::endl;
        // }
        val -= 2 * q; // p[c] reduced to (-2, 2)
    }
};

template<typename T>
struct OptionalFMod2<T, true> {
    static void run(T &/* val */) { /* nop */ }
};

////////////////////////////////////////////////////////////////////////////////
// Symmetry class definitions
////////////////////////////////////////////////////////////////////////////////
// Base unit is a full period cell: [-1, 1]^3
template<typename TOL>
struct TriplyPeriodic : SymmetryCRTP<TriplyPeriodic<TOL>> {
    typedef TOL Tolerance;
    // Disambiguate CRTP instances
    typedef SymmetryCRTP<TriplyPeriodic<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    static constexpr double tolerance = double(TOL::num) / double(TOL::den);
    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() {
        return BBox<Point3<Real>>(Point3<Real>(-1, -1, -1),
                                  Point3<Real>(1, 1, 1));
    }

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        for (size_t c = 0; c < 3; ++c) {
            // Quickly reduce to [-2, 2] for plain scalar types to accelerate
            // far point sampling.
            // This integer conversion-based reduction will not work with
            // autodiff types, but the subsequent brute-force reduction to
            // [-1, 1] should be fast since we should never be querying outside
            // the base cell when using autodiff to compute shape
            // velocities/normals).
            OptionalFMod2<Real>::run(p[c]);

            // Now reduce to [-1, 1] with tolerance...
            while (p[c] >  1.0 + tolerance) p[c] -= 2.0;
            while (p[c] < -1.0 - tolerance) p[c] += 2.0;
        }
        return p;
    }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        return (isNegative<TOL>(std::abs(p[0]) - 1.0)) &&
               (isNegative<TOL>(std::abs(p[1]) - 1.0)) &&
               (isNegative<TOL>(std::abs(p[2]) - 1.0));
    }

    // Note that the group of translational symmetries is infinite, but for our
    // purposes (determining incident edges from neighboring cells), the
    // isometries that take us to adjacent cells are sufficient
    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group(1, Isometry()); // Start with identity element
        // Add in translational symmetries getting us to the 26 adjacent cells
        for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
        for (int z = -1; z <= 1; ++z) {
            if ((x == 0) && (y == 0) && (z == 0)) continue;
            // Note that the base cell is size 2 ([-1, 1]^3)
            group.emplace_back(Isometry::translation(2.0 * x, 2.0 * y, 2.0 * z));
        }
        return group;
    }
};

// Base unit is the positive octant: [0, 1]^3
// Symmetry group D_2h x Translations
template<typename TOL>
struct Orthotropic : public TriplyPeriodic<TOL>, SymmetryCRTP<Orthotropic<TOL>> {
    typedef TOL Tolerance;
    // Disambiguate CRTP instances
    typedef SymmetryCRTP<Orthotropic<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    using TriplyPeriodic<TOL>::tolerance;

    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() {
        return BBox<Point3<Real>>(Point3<Real>(0, 0, 0),
                                  Point3<Real>(1, 1, 1));
    }

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        p = TriplyPeriodic<TOL>::mapToBaseUnit(p);
        for (size_t c = 0; c < 3; ++c)
            if (p[c] < 0) p[c] = -p[c]; // std::abs is problematic for autodiff
        return p;
    }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        return TriplyPeriodic<TOL>::inBaseUnit(p) &&
                isPositive<TOL>(p[0]) && isPositive<TOL>(p[1]) &&
                isPositive<TOL>(p[2]);
    }

    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group;
        std::vector<Isometry> parentGroup = TriplyPeriodic<TOL>::symmetryGroup();
        for (const Isometry &p : parentGroup) {
            group.push_back(p); // Identity reflection

            // Single axis
            group.push_back(p.compose(Isometry::reflection(Axis::X)));
            group.push_back(p.compose(Isometry::reflection(Axis::Y)));
            group.push_back(p.compose(Isometry::reflection(Axis::Z)));

            // Double axis
            group.push_back(p.compose(Isometry::reflection(Axis::X)).compose(Isometry::reflection(Axis::Y)));
            group.push_back(p.compose(Isometry::reflection(Axis::X)).compose(Isometry::reflection(Axis::Z)));
            group.push_back(p.compose(Isometry::reflection(Axis::Y)).compose(Isometry::reflection(Axis::Z)));

            // Triple axis
            group.push_back(p.compose(Isometry::reflection(Axis::X)).compose(Isometry::reflection(Axis::Y)).compose(Isometry::reflection(Axis::Z)));
        }
        
        return group;
    }
};

// Base unit is one of the 6 tetrahedra in the positive octant:
// the tetrahedron with corners (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
// However, we still mesh the positive octant since meshing within a tetrahedron
// is harder.
// Symmetry group Oh x Translations
template<typename TOL>
struct Cubic : public Orthotropic<TOL>, SymmetryCRTP<Cubic<TOL>> {
    typedef TOL Tolerance;
    using Orthotropic<TOL>::representativeMeshCell;
    using TriplyPeriodic<TOL>::tolerance;

    // Disambiguate CRTP instances
    typedef SymmetryCRTP<Cubic<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        p = Orthotropic<TOL>::mapToBaseUnit(p);
        // Sort components descending: an optimal algorithm would still need 3
        // comparisons in the worst case (though 2 in the best), so this bubble
        // sort isn't too bad.
        if (p[0] < p[1]) std::swap(p[0], p[1]);
        if (p[1] < p[2]) std::swap(p[1], p[2]);
        if (p[0] < p[1]) std::swap(p[0], p[1]);
        return p;
    }

    // p is in the canonical base tetrahedron if it's in the positive unit cube
    // and its components are in non-ascending order.
    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        return Orthotropic<TOL>::inBaseUnit(p) &&
                (p[0] + tolerance >= p[1]) &&
                (p[1] + tolerance >= p[2]);
    }

    // Octahedral group symmetries are the reflections of the Orthotropic class'
    // group combined with all 6 axis permutations
    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group;
        std::vector<Isometry> parentGroup = Orthotropic<TOL>::symmetryGroup();
        for (const Isometry &p : parentGroup) {
            group.push_back(p); // X Y Z

            group.push_back(p.compose(Isometry::permutation(Axis::X, Axis::Y))); // Y X Z
            group.push_back(p.compose(Isometry::permutation(Axis::X, Axis::Z))); // Z Y X
            group.push_back(p.compose(Isometry::permutation(Axis::Y, Axis::Z))); // X Z Y

            group.push_back(p.compose(Isometry::permutation(Axis::Y, Axis::Z)).compose(Isometry::permutation(Axis::X, Axis::Y))); // Y Z X
            group.push_back(p.compose(Isometry::permutation(Axis::X, Axis::Y)).compose(Isometry::permutation(Axis::Y, Axis::Z))); // Z X Y
        }
        return group;
    }
};

} // end of namespace Symmetry

#endif /* end of include guard: SYMMETRY_HH */

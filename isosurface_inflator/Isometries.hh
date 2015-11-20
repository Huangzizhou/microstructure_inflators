////////////////////////////////////////////////////////////////////////////////
// Isometries.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Represents an isometry, intended to be used as an element of a symmetry
//      group. These can also be used to map points.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/26/2015 18:01:53
////////////////////////////////////////////////////////////////////////////////
#ifndef ISOMETRIES_HH
#define ISOMETRIES_HH
#include <memory>
#include <iostream>
#include <cassert>
#include <string>

namespace Symmetry {
    enum class Axis : unsigned int { X = 0, Y = 1, Z = 2 };
    inline std::string axisName(Axis a) {
        switch(a) {
            case Axis::X: return "X";
            case Axis::Y: return "Y";
            case Axis::Z: return "Z";
        }
        return "Unknown";
    }
}

// Elements of a symmetry group
struct Isometry {
    struct Operation;

    // Default constructor: identity isometry
    Isometry() { }
    Isometry(std::shared_ptr<const Operation> op) : operations(1, op) { }

    bool isIdentity() const { return operations.empty(); }

    // Compose "op" on the right (op will be performed before operations[]).
    Isometry compose(std::shared_ptr<const Operation> op) const {
        Isometry result(*this);
        result.operations.push_back(op);
        return result;
    }

    // Performed in right-to-left: operations[n - 1], ..., operations[0]
    template<typename Real>
    Point3<Real> apply(Point3<Real> p) const {
        for (size_t i = operations.size(); i > 0; --i) {
            operations[i - 1]->apply(p);
        }
        return p;
    }

    void print(std::ostream &os) const {
        bool first = true;
        if (operations.size() == 0) os << "identity";
        for (auto &op : operations) {
            if (!first) os << " * ";
            op->print(os);
            first = false;
        }
    }

    struct Operation {
        template<typename Real>
        void apply(Point3<Real> &p) const { return applyOperation(this, p); }
        virtual void print(std::ostream &os) const = 0;
        virtual ~Operation() { }
    };

    struct Translation : public Operation {
        Translation(double x, double y, double z) : t(x, y, z) { }
        virtual ~Translation() { }
        virtual void print(std::ostream &os) const { os << "t(" << t[0] << ", " << t[1] << ", " << t[2] << ")"; }
        Vector3<double> t;
    };

    struct Reflection : public Operation {
        Reflection(Symmetry::Axis a) : a(a) { }
        virtual ~Reflection() { }
        virtual void print(std::ostream &os) const { os << "r(" << axisName(a) << ")"; }
        Symmetry::Axis a;
    };

    struct Permutation : public Operation {
        Permutation(Symmetry::Axis a1, Symmetry::Axis a2) : a1(a1), a2(a2) { }
        virtual ~Permutation() { }
        virtual void print(std::ostream &os) const { os << "p(" << axisName(a1) << ", " << axisName(a2) << ")"; }
        Symmetry::Axis a1, a2;
    };

    // Unfortunately virtual template methods are not allowed, so we're forced
    // to implement polymorphism with a dynamic cast in applyOperation()
    template<typename Real>
    static void applyOperation(const Operation *op, Point3<Real> &p) {
        if (const Translation *trans = dynamic_cast<const Translation *>(op)) {
            for (size_t i = 0; i < 3; ++i) p[i] += trans->t[i];
        }
        else if (const Reflection *refl = dynamic_cast<const Reflection *>(op)) {
            p[static_cast<unsigned int>(refl->a)] *= -1;
        }
        else if (const Permutation *perm = dynamic_cast<const Permutation *>(op)) {
            std::swap(p[static_cast<unsigned int>(perm->a1)],
                      p[static_cast<unsigned int>(perm->a2)]);
        }
        else { assert("Unrecognized operation!"); } // Impossible
    }

    static std::shared_ptr<Translation> translation(double x, double y, double z)         { return std::make_shared<Translation>(x, y, z); }
    static std::shared_ptr< Reflection>  reflection(Symmetry::Axis a)                     { return std::make_shared< Reflection>(a); }
    static std::shared_ptr<Permutation> permutation(Symmetry::Axis a1, Symmetry::Axis a2) { return std::make_shared<Permutation>(a1, a2); }

    std::vector<std::shared_ptr<const Operation>> operations;
};

inline std::ostream &operator<<(std::ostream &os, const Isometry &i) {
    i.print(os);
    return os;
}

#endif /* end of include guard: ISOMETRIES_HH */

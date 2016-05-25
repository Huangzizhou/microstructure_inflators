////////////////////////////////////////////////////////////////////////////////
// IterateBase.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Non-templated, dynamic interface to the pattern optimization iterate.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  05/02/2016 19:08:21
////////////////////////////////////////////////////////////////////////////////
#ifndef ITERATEBASE_HH
#define ITERATEBASE_HH

#include <vector>
#include <memory>
#include <stdexcept>

#include <Fields.hh>
#include "EvaluatedObjectiveTerm.hh"

namespace PatternOptimization {

struct IterateBase {
    using EOTPtr = std::unique_ptr<EvaluatedObjectiveTerm>;

    using SField = ScalarField<Real>;
    // parametricOptimization: whether we're running a lower-dimensional
    // parametric optimization where gradients must be computed by taking inner
    // products with the shape velocity fields.
    // 
    // For free boundary shape optimization problems, this inner product
    // approach to gradients is intractable and unnecessary.
    IterateBase(bool parametricOptimization) : m_parametricOptimization(parametricOptimization) { }
    const EvaluatedObjectiveTerm &evaluatedObjectiveTerm(const std::string &name) const {
        m_assertParametric();
        for (const auto &term : m_evaluatedObjectiveTerms)
            if (term->name == name) return *term;
        throw std::runtime_error("Objective term not found: " + name);
    }

    const EvaluatedObjectiveTerm &evaluatedObjectiveTerm(size_t i) const {
        m_assertParametric();
        return *m_evaluatedObjectiveTerms.at(i);
    }

    const std::vector<EOTPtr> &evaluatedObjectiveTerms() const {
        m_assertParametric();
        return m_evaluatedObjectiveTerms;
    }

    virtual size_t numObjectiveTerms() const { return m_evaluatedObjectiveTerms.size(); }

    ////////////////////////////////////////////////////////////////////////////
    // Full objective parametric gradient and steepest descent: linear
    // combination of evaluated objective terms' gradients.
    ////////////////////////////////////////////////////////////////////////////
    SField gradp() const {
        SField full;
        for (const auto &term : m_evaluatedObjectiveTerms) {
            if (full.domainSize() == 0) full  = term->gradContribution();
            else                        full += term->gradContribution();
        }
        return full;
    }

    SField steepestDescentParam() const {
        m_assertParametric();
        // TODO: construct parameter "mass matrix"...
        // Normalize for unit M-norm
        throw std::runtime_error("Unimplemented.");
    }

    bool paramsDiffer(size_t nParams, const Real *params) const {
        assert(nParams = m_params.size());
        for (size_t i = 0; i < nParams; ++i)
            if (m_params[i] != params[i])
                return true;
        return false;
    }

    const std::vector<Real> &params() const { return m_params; }

    // Implemented by subclass
    virtual Real evaluate() const = 0;
    virtual void writeMeshAndFields(const std::string &path) const = 0;
    virtual void writeDescription(std::ostream &os) const = 0;
    virtual std::vector<std::string> objectiveTermNames() const = 0;

    bool isParametric() const { return m_parametricOptimization; }

protected:
    // Current params
    std::vector<Real> m_params;

    // Whether we're running a lower-dimensional parametric optimization where
    // gradients must be computed by taking inner products with the shape
    // velocity fields. This also determines whether gradients are printed by
    // the writeDescription method.
    bool m_parametricOptimization;

    void m_assertParametric() const {
        if (!m_parametricOptimization)
            throw std::runtime_error("This operation is only supported for parametric optimization.");
    }

    // Filled out by subclass
    std::vector<EOTPtr> m_evaluatedObjectiveTerms;
};

}

#endif /* end of include guard: ITERATEBASE_HH */
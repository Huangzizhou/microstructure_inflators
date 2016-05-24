////////////////////////////////////////////////////////////////////////////////
// PatternOptimizationIterate.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Encapsulates the state of a pattern optimization iterate and provides
//      objective/gradient/etc.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/26/2014 19:04:20
////////////////////////////////////////////////////////////////////////////////
#ifndef PATTERNOPTIMIZATIONITERATE_HH
#define PATTERNOPTIMIZATIONITERATE_HH

#include <EdgeFields.hh>
#include <MSHFieldWriter.hh>

#include <iostream>
#include <cstdio>
#include <cassert>
#include <memory>
#include <iostream>
#include <iomanip>
#include <memory>
#include <tuple>

#include <MeshIO.hh>
#include <TriMesh.hh>
#include <filters/subdivide.hh>

#include "PatternOptimizationConfig.hh"
#include "ObjectiveTerm.hh"
#include "EvaluatedObjectiveTerm.hh"
#include "SDConversions.hh"

#include "IterateBase.hh"

#include <Future.hh>

namespace PatternOptimization {

template<class _Sim>
struct Iterate : public IterateBase {
    using   VField = typename _Sim::VField;
    using   SField = ScalarField<Real>;
    using _ETensor = typename _Sim::ETensor;
    static constexpr size_t _N = _Sim::N;

    using ObjectiveTermPtr = std::unique_ptr<ObjectiveTerm<_N>>;
    using ObjectiveTermMap = std::map<std::string, ObjectiveTermPtr>;

    using IterateBase::m_params;

    Iterate(Inflator<_N> &inflator, size_t nParams, const double *params)
        : IterateBase(inflator.isParametric())
    {
        m_params.resize(nParams);
        for (size_t i = 0; i < nParams; ++i)
            m_params[i] = params[i];
        m_printable = inflator.isPrintable(m_params);

        // std::cout << "Inflating" << std::endl;
        BENCHMARK_START_TIMER_SECTION("Inflate");
        try {
            inflator.inflate(m_params);
        }
        catch (...) {
            // Hack to correct timer behavior--should probably use RAII
            BENCHMARK_STOP_TIMER_SECTION("Inflate");
            throw;
        }
        BENCHMARK_STOP_TIMER_SECTION("Inflate");
        // std::cout << "Inflated" << std::endl;

        // std::cout << "Checking geometry" << std::endl;
        if ((inflator.elements().size() == 0) || (inflator.vertices().size() == 0)) {
            std::cerr << std::setprecision(20);
            std::cerr << "Exception while inflating parameters" << std::endl;
            for (size_t i = 0; i < m_params.size(); ++i) std::cerr << m_params[i] << "\t";
            std::cerr << std::endl;
            throw std::runtime_error("Empty inflated geometry. Elements: "
                    + std::to_string(inflator.elements().size()) + ", Vertices: "
                    + std::to_string(inflator.vertices().size()));
        }
        // std::cout << std::endl;


        // std::cout << "Building Simulator" << std::endl;
        BENCHMARK_START_TIMER_SECTION("Eval");

        {
            std::unique_ptr<std::vector<MeshIO::IOVertex >> verts;
            std::unique_ptr<std::vector<MeshIO::IOElement>> elems;
            m_sim = Future::make_unique<_Sim>(inflator.elements(),
                                              inflator.vertices());
        }
        // std::cout << "Done" << std::endl;
        // std::cout << "Homogenizing" << std::endl;

        try {
            PeriodicHomogenization::solveCellProblems(w_ij, *m_sim, 1e-9);
        }
        catch(std::exception &e) {
            std::cerr << "Cell problem solve failed: " << e.what() << std::endl;
            MeshIO::save("debug.msh", mesh());
            std::cerr << "Wrote geometry to 'debug.msh'" << std::endl;
            std::cerr << std::setprecision(19) << std::endl;
            std::cerr << "params:";
            for (size_t i = 0; i < m_params.size(); ++i) std::cerr << "\t" << m_params[i];
            std::cerr << std::endl;
            exit(-1);
        }

        C = PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(w_ij, *m_sim);
        S = C.inverse();

        // std::cout << "Done" << std::endl;

        BENCHMARK_STOP_TIMER_SECTION("Eval");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Convenience methods for evaluating the full objective and its shape
    // derivatives/gradients/steepest descent directions
    ////////////////////////////////////////////////////////////////////////////////
    // Evaluate full objective
    Real evaluate() const {
        Real full = 0;
        if (isParametric()) {
            for (const auto &term : m_evaluatedObjectiveTerms)
                full += term->contribution();
        }
        else {
            for (const auto &term : m_objectiveTerms)
                full += term.second->normalizedWeight() * term.second->evaluate();
        }
        return full;
    }

    // Evaluate normalized sub-objective
    Real evaluateNormalized(const std::string &name) const {
        return m_objectiveTerms.at(name)->evaluateNormalized();
    }

    // Evaluate full objective's differential one-form acting on boundary
    // velocity vector fields
    VField differential() const {
        VField full;
        for (const auto &term : m_objectiveTerms) {
            VField contrib = term.second->diff_bdry();
            contrib *= term.second->normalizedWeight();
            if (full.domainSize() == 0) full  = contrib;
            else                        full += contrib;
        }
        return full;
    }

    // Steepest descent boundary velocity field for the full objective
    VField steepestDescent() const {
        return SDConversions::descent_from_diff_bdry(differential(), m_sim);
    }

    // Steepest descent parameter vector in terms of object's L2 metric. 
    SField steepestDescentParam(const std::vector<VField> &/*bdrySVels*/) const {
        throw std::runtime_error("steepestDescentParam unimplemented");
    }

    // Partial derivatives of full objective wrt each pattern parameter with
    // specified boundary velocities.
    using IterateBase::gradp;   // Prevent hiding over overload
    SField gradp(const std::vector<VField> &bdrySVels) const {
        auto diff = differential();
        SField g(bdrySVels.size());
        for (size_t p = 0; p < bdrySVels.size(); ++p) {
            assert(bdrySVels[p].domainSize() == diff.domainSize());
            g[p] = diff.innerProduct(bdrySVels[p]);
        }
        return g;
    }

    const ObjectiveTerm<_N> &objectiveTerm(const std::string &name) const { return *m_objectiveTerms.at(name); }
          ObjectiveTerm<_N> &objectiveTerm(const std::string &name)       { return *m_objectiveTerms.at(name); }
    const ObjectiveTermMap  &objectiveTerms() const { return m_objectiveTerms; }

    virtual size_t numObjectiveTerms() const { return m_objectiveTerms.size(); }

    virtual std::vector<std::string> objectiveTermNames() const {
        std::vector<std::string> names;
        for (auto &term : m_objectiveTerms) names.push_back(term.first);
        return names;
    }

    void addObjectiveTerm(const std::string &name, std::unique_ptr<ObjectiveTerm<_N>> t) {
        if (m_objectiveTerms.count(name))
            throw std::runtime_error("Objective term '" + name + "'already added.");
        m_objectiveTerms.emplace(name, std::move(t));
    }

    // Estimate all (evaluated) sub-objectives at an offset point: used for
    // approximately evaluating uninflatable points.
    //
    // WARNING: this affects the result of the EvaluatedObjectiveTerms' value() and
    // this->evaluate() for parametric optimization iterates, but not the
    // individual objective terms' evaluate().
    void estimatePoint(size_t nParams, const double *params) {
        this->m_assertParametric();
        assert(nParams == m_params.size());
        SField delta(nParams);
        for (size_t p = 0; p < nParams; ++p)
            delta[p] = params[p] - m_params[p];
        for (auto &term : m_evaluatedObjectiveTerms)
            term->setEstimateWithDeltaParams(delta);

        std::cerr << "WARNING, USING APPROXIMATE OBJECTIVES/GRADIENTS AT DIST:";
        for (size_t p = 0; p < delta.domainSize(); ++p) std::cerr << "\t" << delta[p];
        std::cerr << std::endl;
    }

    void disableEstimation() {
        for (auto &term : m_evaluatedObjectiveTerms)
            term->disableEstimation();
    }

    const _Sim &simulator() const { return *m_sim; }
          _Sim &simulator()       { return *m_sim; }

          typename _Sim::Mesh &mesh()       { return m_sim->mesh(); }
    const typename _Sim::Mesh &mesh() const { return m_sim->mesh(); }

    const _ETensor &elasticityTensor() const { return C; }
    const _ETensor &complianceTensor() const { return S; }
    const std::vector<VField> &fluctuationDisplacements() const { return w_ij; }

    void evaluateObjectiveTerms(const Inflator<_N> &inflator) {
        m_evaluatedObjectiveTerms.clear();
        m_evaluatedObjectiveTerms.reserve(m_objectiveTerms.size());
        for (auto &term : m_objectiveTerms) {
            std::unique_ptr<EvaluatedObjectiveTerm> eterm;
            if (auto nllsTerm = dynamic_cast<const NLLSObjectiveTerm<_N> *>(term.second.get())) {
                auto enterm = Future::make_unique<EvaluatedObjectiveTermNLLS>();
                // Fill out NLLS info
                enterm->residualComponents = nllsTerm->residual();
                enterm->jacobianComponents = nllsTerm->jacobian(inflator.shapeVelocities(mesh()));
                eterm = std::move(enterm);
            }
            else {
                eterm = Future::make_unique<EvaluatedObjectiveTerm>();
            }
            assert(eterm);
            eterm->name = term.first;
            eterm->setValue(term.second->evaluate());

            eterm->normalization = term.second->normalization();
            eterm->weight        = term.second->weight();
            eterm->gradp         = term.second->gradp(inflator.shapeVelocities(mesh()));

            IterateBase::m_evaluatedObjectiveTerms.push_back(std::move(eterm));
        }
    }

    void writeMeshAndFields(const std::string &path) const {
        MSHFieldWriter writer(path, m_sim->mesh());
        for (auto &term : m_objectiveTerms)
            term.second->writeFields(writer);
    }

    void writeDescription(std::ostream &os) const {
        if (isParametric()) {
            os << "p:";
            for (size_t i = 0; i < m_params.size(); ++i)
                os << "\t" << m_params[i];
            os << std::endl;
        }

        os << "moduli:\t";
        C.printOrthotropic(os);
        os << "anisotropy:\t" << C.anisotropy() << std::endl;
        os << "printable:\t" << m_printable << std::endl;

        // TODO: non-parametric gradient norm info:
        // M_norm(steepestDescent), since steepest descent is the Riesz representative of
        // the differential, and we want it's norm. This ends up being
        // sqrt(g^T M^-1 g) where g is the differential

        for (auto &term : m_objectiveTerms)
            term.second->writeDescription(os, term.first);
        if (isParametric()) {
            // Evaluated objective terms know the gradient information
            for (auto &eterm : this->m_evaluatedObjectiveTerms)
                eterm->writeGradientDescription(os);
        }

        if (numObjectiveTerms() > 1) {
            os << "JFull:\t" << this->evaluate() << std::endl;
            if (isParametric()) {
                SField gp = IterateBase::gradp();
                os << "grad_p JFull:\t";
                gp.print(os, "", "", "", "\t");
                os << std::endl << "||grad_p JFull||:\t" << gp.norm() << std::endl;
            }
        }
    }

protected:
    std::unique_ptr<_Sim> m_sim;
    _ETensor C, S ;
    bool m_printable;

    bool m_parametricOptimization;

    ObjectiveTermMap m_objectiveTerms; 

    // Fluctuation displacements
    std::vector<VField> w_ij;
};

}

#endif /* end of include guard: PATTERNOPTIMIZATIONITERATE_HH */

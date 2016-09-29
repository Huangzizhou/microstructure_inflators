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
#include "BaseCellOperations.hh"

#include "IterateBase.hh"
#include "Inflator.hh"

#include <Future.hh>

namespace PatternOptimization {

template<class _Sim>
struct Iterate : public IterateBase {
    using    OForm = ScalarOneForm<_Sim::N>;
    using   VField = typename _Sim::VField;
    using   SField = ScalarField<Real>;
    using _ETensor = typename _Sim::ETensor;
    static constexpr size_t _N = _Sim::N;

    using ObjectiveTermPtr = std::unique_ptr<ObjectiveTerm<_N>>;
    using ObjectiveTermMap = std::map<std::string, ObjectiveTermPtr>;

    using IterateBase::m_params;

    // Ortho base cell only? TODO: determine from inflator
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

        // MeshIO::save("debug.msh", inflator.vertices(), inflator.elements());
        // if (isParametric()) {
        //     std::cerr << std::setprecision(19) << std::endl;
        //     std::cerr << "params:";
        //     for (size_t i = 0; i < m_params.size(); ++i) std::cerr << "\t" << m_params[i];
        // }

        // std::cout << "Building Simulator" << std::endl;
        BENCHMARK_START_TIMER_SECTION("Eval");

        m_sim = Future::make_unique<_Sim>(inflator.elements(),
                                          inflator.vertices());
        // std::cout << "Done" << std::endl;
        // std::cout << "Homogenizing" << std::endl;

        try {
            m_baseCellOps = BaseCellOperations<_Sim>::construct(inflator.baseCellType(), *m_sim);
        }
        catch(std::exception &e) {
            std::cerr << "Cell problem solve failed: " << e.what() << std::endl;
            MeshIO::save("debug.msh", mesh());
            std::cerr << "Wrote geometry to 'debug.msh'" << std::endl;
            std::cerr << std::setprecision(19) << std::endl;
            if (isParametric()) {
                std::cerr << "params:";
                for (size_t i = 0; i < m_params.size(); ++i) std::cerr << "\t" << m_params[i];
            }
            std::cerr << std::endl;
            exit(-1);
        }

        C = m_baseCellOps->homogenizedElasticityTensor();
        S = C.inverse();

        // std::cout << "Done" << std::endl;

        BENCHMARK_STOP_TIMER_SECTION("Eval");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Convenience methods for evaluating the full objective and its shape
    // derivatives/gradients/steepest descent directions
    ////////////////////////////////////////////////////////////////////////////////
    // Evaluate full objective
    virtual Real evaluate() const override {
        Real full = 0;
        for (const auto &term : m_evaluatedObjectiveTerms)
            full += term->contribution();
        return full;
    }

    // Evaluate normalized sub-objective
    virtual Real evaluateNormalized(const std::string &name) const override {
        return m_objectiveTerms.at(name)->evaluateNormalized();
    }

    // Evaluate full objective's differential one-form acting on boundary
    // velocity vector fields
    OForm differential() const {
        OForm full;
        for (const auto &term : m_objectiveTerms) {
            OForm contrib = term.second->diff_bdry();
            contrib *= term.second->normalizedWeight();
            if (full.domainSize() == 0) full  = contrib;
            else                        full += contrib;
        }
        return full;
    }

    // Partial derivatives of full objective wrt each pattern parameter with
    // specified boundary velocities.
    using IterateBase::gradp;   // Prevent hiding instead of overloading
    SField gradp(const std::vector<VField> &bdrySVels) const {
        OForm diff = differential();
        SField g(bdrySVels.size());
        for (size_t p = 0; p < bdrySVels.size(); ++p)
            g[p] = diff[bdrySVels[p]];
        return g;
    }

    const ObjectiveTerm<_N> &objectiveTerm(const std::string &name) const { return *m_objectiveTerms.at(name); }
          ObjectiveTerm<_N> &objectiveTerm(const std::string &name)       { return *m_objectiveTerms.at(name); }
    const ObjectiveTermMap  &objectiveTerms() const { return m_objectiveTerms; }

    virtual size_t numObjectiveTerms() const override { return m_objectiveTerms.size(); }

    virtual std::vector<std::string> objectiveTermNames() const override {
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
        assert(nParams == m_params.size());
        SField delta(nParams);
        for (size_t p = 0; p < nParams; ++p)
            delta[p] = params[p] - m_params[p];
        for (auto &term : m_evaluatedObjectiveTerms)
            term->setEstimateWithDeltaParams(delta);

        std::cerr << "WARNING, USING APPROXIMATE OBJECTIVES/GRADIENTS AT DIST:";
        if (isParametric())
            for (size_t p = 0; p < delta.domainSize(); ++p) std::cerr << "\t" << delta[p];
        else delta.norm();
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

          BaseCellOperations<_Sim> &baseCellOps()       { return *m_baseCellOps; }
    const BaseCellOperations<_Sim> &baseCellOps() const { return *m_baseCellOps; }

    const _ETensor &elasticityTensor() const { return C; }
    const _ETensor &complianceTensor() const { return S; }
    const std::vector<VField> &fluctuationDisplacements() const { return m_baseCellOps->fluctuationDisplacements(); }

    // Must be called *after* objective terms are added. (cannot be from the
    // constructor...)
    void evaluateObjectiveTerms(const Inflator<_N> &inflator) {
        m_evaluatedObjectiveTerms.clear();
        m_evaluatedObjectiveTerms.reserve(m_objectiveTerms.size());
        for (auto &term : m_objectiveTerms) {
            std::unique_ptr<EvaluatedObjectiveTerm> eterm;
            // Conditionally construct EvaluatedObjectiveTerm{,NLLS}
            {
                auto nllsTerm = dynamic_cast<const NLLSObjectiveTerm<_N> *>(term.second.get());
                // NLLS optimization only works for parametric inflators--use
                // general nonlinear optimization otherwise
                if (inflator.isParametric() && nllsTerm) {
                    auto enterm = Future::make_unique<EvaluatedObjectiveTermNLLS>();
                    // Fill out NLLS info
                    enterm->residualComponents = nllsTerm->residual();
                    enterm->jacobianComponents = nllsTerm->jacobian(inflator.shapeVelocities(mesh()));
                    eterm = std::move(enterm);
                }
                else {
                    eterm = Future::make_unique<EvaluatedObjectiveTerm>();
                }
            }
            assert(eterm);
            eterm->name = term.first;
            eterm->setValue(term.second->evaluate());

            eterm->normalization = term.second->normalization();
            eterm->weight        = term.second->weight();
            if (inflator.isParametric()) {
                eterm->gradp = term.second->gradp(inflator.shapeVelocities(mesh()));
                // TODO: construct parameter "mass matrix"...
                // Normalize for unit M-norm
                eterm->descentp = eterm->gradp;
                eterm->descentp *= -1.0;
            }
            else {
                // TODO: use baseCellOps for this... (not needed for parametric
                // optimization).
                std::cerr << "Computing grad for " << term.first << std::endl;
                eterm->gradp = inflator.paramsFromBoundaryVField(term.second->differential().asVectorField());
                std::cerr << "Computing descent for " << term.first << std::endl;
                eterm->descentp = inflator.paramsFromBoundaryVField(
                        SDConversions::descent_from_diff_bdry(term.second->differential(), *m_sim));
                std::cerr << "Done with " << term.first << std::endl;
            }

            IterateBase::m_evaluatedObjectiveTerms.push_back(std::move(eterm));
        }
    }

    virtual void writeMeshAndFields(const std::string &path) const override {
        MSHFieldWriter writer(path, m_sim->mesh());
        for (auto &term : m_objectiveTerms)
            term.second->writeFields(writer);
    }

    virtual void writeDescription(std::ostream &os) const override {
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
        // Evaluated objective terms know the gradient information
        for (auto &eterm : this->m_evaluatedObjectiveTerms)
            eterm->writeGradientDescription(os, isParametric());

        if (numObjectiveTerms() > 1) {
            os << "JFull:\t" << this->evaluate() << std::endl;
            SField gp = IterateBase::gradp();
            if (isParametric()) {
                os << "grad_p JFull:\t";
                gp.print(os, "", "", "", "\t");
                os << std::endl;
            }
            os << "||grad_p JFull||:\t" << gp.norm() << std::endl;
        }
    }

    virtual ~Iterate() { }

protected:
    std::unique_ptr<_Sim> m_sim;
    std::unique_ptr<BaseCellOperations<_Sim>> m_baseCellOps;
    _ETensor C, S ;
    bool m_printable;

    ObjectiveTermMap m_objectiveTerms; 
};

}

#endif /* end of include guard: PATTERNOPTIMIZATIONITERATE_HH */

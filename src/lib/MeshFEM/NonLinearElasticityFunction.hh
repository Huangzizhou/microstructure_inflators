//
// Created by Davi Colli Tozoni on 7/19/18.
//

#ifndef NONLINEARELASTICITYFUNCTION_H
#define NONLINEARELASTICITYFUNCTION_H

#include <MeshFEM/SparseMatrices.hh>

template<typename Real>
class SystemTransformations {
public:
    SystemTransformations() {}

    SystemTransformations(size_t originalSize, const std::vector<size_t> &fixedVars, const std::vector<Real>  &fixedVarValues) :
            m_fixedVars(fixedVars), m_fixedVarValues(fixedVarValues) {

        assert(fixedVars.size() == fixedVarValues.size());

        // Now, it is possible to run without any fixed vars.
        //if (fixedVars.size() == 0) return;

        m_reducedVarForVar.resize(originalSize);
        // Identity mapping of variables to reduced variables.
        for (size_t i = 0; i < originalSize; ++i)
            m_reducedVarForVar[i] = i;
        m_fixedVarValues.clear();

        // replacementIndex tracks what the current reduced variable indices are
        // remapped to. Initially it is used to flag (reduced) variables for
        // elimination (with -1), but afterward the full array is filled in.
        std::vector<int> replacementIndex(originalSize, 0);

        // The value to which each (reduced) variable will be fixed, or zero if
        // the variable will not be fixed. Needed for efficiently computing RHS
        // contribution of fixedVarValues
         m_originalIndexToFixedValues = std::vector<Real>(originalSize, 0.0);
        for (size_t i = 0; i < fixedVars.size(); ++i) {
            int rv = m_reducedVarForVar[fixedVars[i]];
            if (rv < 0) continue;
            m_originalIndexToFixedValues[rv] = fixedVarValues[i];
        }

        // Mark fixed variables for elimination and store their values
        // m_fixedVarValues for post-solve recovery.
        // Also move fixedVarValues[i] over to m_reducedVarForVar
        for (size_t i = 0; i < fixedVars.size(); ++i) {
            size_t toFix = fixedVars[i];
            assert(toFix < m_reducedVarForVar.size());

            // Get the current reduced index of the variable.
            int curr = m_reducedVarForVar[toFix];
            if (curr < 0) throw std::runtime_error("Variable already fixed.");
            assert(size_t(curr) < replacementIndex.size());

            replacementIndex[curr] = -1;
            m_reducedVarForVar[toFix] = -1 - int(m_fixedVarValues.size());
            Real val = fixedVarValues[i];
            m_fixedVarValues.push_back(val);
        }

        // Reindex all the current reduced variables.
        size_t newIdx = 0;
        for (size_t i = 0; i < originalSize; ++i) {
            if (replacementIndex[i] >= 0)
                replacementIndex[i] = newIdx++;
        }

        // Apply replacement to m_reducedVarForVar.
        for (size_t i = 0; i < originalSize; ++i) {
            int curr = m_reducedVarForVar[i];
            if (curr < 0) continue;
            assert(size_t(curr) < replacementIndex.size());
            m_reducedVarForVar[i] = replacementIndex[curr];
        }
    }

    std::vector<Real> reducedDisplacementToOriginalVector(std::vector<Real> &uReduced) const {
        std::vector<Real> u(m_reducedVarForVar.size());

        for (size_t v = 0; v < m_reducedVarForVar.size(); v++) {
            int r = m_reducedVarForVar[v];
            if (r < 0) {
                size_t fixedVar = -1 - r;
                assert(fixedVar < m_fixedVarValues.size());
                u[v] = m_fixedVarValues[fixedVar];
            }
            else {
                assert(size_t(r) < uReduced.size());
                u[v] = uReduced[r];
            }
        }

        return u;
    }

    std::vector<Real> reducedToOriginalVector(std::vector<Real> &reduced) const {
        std::vector<Real> result(m_reducedVarForVar.size());

        for (size_t v = 0; v < m_reducedVarForVar.size(); v++) {
            int r = m_reducedVarForVar[v];
            if (r < 0)
                result[v] = 0.0;
            else {
                assert(size_t(r) < reduced.size());
                result[v] = reduced[r];
            }
        }

        return result;
    }

    template<class _Vec>
    std::vector<Real> originalToReducedVector(_Vec &original) const {
        std::vector<Real> reduced(original.size() - m_fixedVarValues.size());

        for (size_t v = 0; v < original.size(); v++) {
            int r = m_reducedVarForVar[v];
            if (r < 0)
                continue;

            reduced[r] = original[v];
        }

        return reduced;
    }

    TripletMatrix<Triplet<Real>> originalToReducedMatrix(TripletMatrix<Triplet<Real>> &original) const {
        TripletMatrix<Triplet<Real>> reduced = original;

        // Remove entries (rows, cols) of original matrix
        auto newEnd = std::remove_if(reduced.nz.begin(), reduced.nz.end(),
                                     [&](const Triplet<Real> &t) -> bool {
                                         return (m_reducedVarForVar[t.i] < 0) ||
                                                (m_reducedVarForVar[t.j] < 0); });
        reduced.nz.erase(newEnd, reduced.nz.end());

        // Apply replacement to A matrix.
        for (Triplet<Real> &t : reduced.nz) {
            t.i = m_reducedVarForVar[t.i];
            t.j = m_reducedVarForVar[t.j];
        }

        reduced.m -= m_fixedVarValues.size();
        reduced.n -= m_fixedVarValues.size();

        return reduced;
    }

    size_t originalSize() const {
        return m_originalIndexToFixedValues.size();
    }


    std::vector<int> m_reducedVarForVar;
    std::vector<size_t> m_fixedVars;
    std::vector<Real> m_fixedVarValues;
    std::vector<Real> m_originalIndexToFixedValues;
};

template<typename Real>
class NonLinearElasticityFunction {
public:
    virtual ~NonLinearElasticityFunction() = default;

    // Non-reduced form: meaning that u corresponds to displacement in all nodes. Also, result is nonlinear term for
    // all nodes of the mesh
    virtual std::vector<Real> evaluate(std::vector<Real> u) const = 0;

    // Non-reduced form: meaning that u corresponds to displacement in all nodes. Also, result is jacobian matrix for
    // all nodes of the mesh, which means matrix has size (nodes mesh x dimension)^2.
    virtual TripletMatrix<Triplet<Real>> jacobian(std::vector<Real> u) const = 0;

    // Non-reduced form: meaning that u corresponds to displacement in all nodes. Also, result is jacobian matrix for
    // all nodes of the mesh, which means matrix has size (nodes mesh x dimension)^2.
    virtual TripletMatrix<Triplet<Real>> approximateJacobian(std::vector<Real> u) const {
        std::vector<Real> perturbedU = u;
        std::vector<Real> negPerturbedU = u;
        Real perturbation = 1e-9;
        TripletMatrix<Triplet<Real>> approximatedJacobian(u.size(), u.size());

        // Evaluate function at u
        std::vector<Real> result = evaluate(u);

        for (unsigned j = 0; j < u.size(); j++) {
            // Compute perturbed u
            perturbedU[j] = u[j] + perturbation;
            negPerturbedU[j] = u[j] - perturbation;

            // Evaluate function on perturbed u
            std::vector<Real> perturbedResult = evaluate(perturbedU);
            std::vector<Real> negPerturbedResult = evaluate(negPerturbedU);

            //std::cout << "result.size()" << result.size() << std::endl;
            for (unsigned i = 0; i < result.size(); i++) {
                Real relativeDifference = (perturbedResult[i] - negPerturbedResult[i]) / (2*perturbation);

                //std::cout << "i: " << i << std::endl;
                //std::cout << "j: " << j << std::endl;
                //std::cout << "result = " << result[i] << std::endl;
                if (abs(relativeDifference) > 1e-20) {
                    //std::cout << "result = " << result[i] << std::endl;
                    //std::cout << "perturbed = " << perturbedResult[i] << std::endl;
                    //std::cout << "Relative difference = " << relativeDifference << std::endl;
                    approximatedJacobian.addNZ(i, j, relativeDifference);
                }
            }

            perturbedU[j] = u[j];
            negPerturbedU[j] = u[j];
        }

        return approximatedJacobian;
    }

    bool compareJacobians(TripletMatrix<Triplet<Real>> &jacobian, TripletMatrix<Triplet<Real>> &approximatedJacobian) const {
        TripletMatrix<Triplet<Real>> diff(jacobian.m, jacobian.n);

        // Compute matrix which represents the difference between the two jacobians
        for (auto t : jacobian.nz) {
            diff.addNZ(t.i, t.j, t.v);
        }
        // Add each element of negative approximated Jacobian to actual jacobian matrix. And sum repeated.
        // Result should be empty matrix
        for (auto t : approximatedJacobian.nz) {
            diff.addNZ(t.i, t.j, -t.v);
        }
        diff.sumRepeated();

        // Verify that all entries are super small or zero
        int errorCount = 0;
        for (auto t : diff.nz) {
            if (abs(t.v) > 1e-8) {
                Triplet<Real> originalTriplet = jacobian.getTriplet(t.i, t.j);
                std::cerr << "[Warning!!] Element (" << t.i << ", " << t.j << ") differs in jacobian by " << t.v
                          << " (original value: " << originalTriplet.v << ")" << std::endl;
                errorCount++;
            }
        }

        //std::cout << "Number of nz elements in Jacobian: " << jacobian.nz.size() << std::endl;
        //std::cout << "Number of nz elements in approximate Jacobian: " << approximatedJacobian.nz.size() << std::endl;
        if (errorCount > 0)
        std::cerr << "Error count: " << errorCount << std::endl;

        bool result = errorCount > 0 ? false : true;

        return result;
    }

    // Size. Number of vector positions
    virtual size_t size() const = 0;
};

template<typename Real>
class ReducedNonLinearElasticityFunction : public NonLinearElasticityFunction<Real> {
public:

    ReducedNonLinearElasticityFunction(std::shared_ptr<NonLinearElasticityFunction<Real>> originalFunction) :
            m_function(originalFunction){
        m_size = 0;
    }

    ReducedNonLinearElasticityFunction(NonLinearElasticityFunction<Real> &originalFunction,
                                       std::vector<size_t> &fixedVars, const std::vector<Real>  &fixedVarValues) :
            m_function(originalFunction), m_systemTransformations(originalFunction.size(), fixedVars, fixedVarValues) {
        m_size = originalFunction.size() - fixedVars.size();
    }

    void fixVariables(const std::vector<size_t> &fixedVars, const std::vector<Real>  &fixedVarValues) {
        m_systemTransformations = SystemTransformations<Real>(m_function->size(), fixedVars, fixedVarValues);
    }

    // Reduced form: meaning that u corresponds to displacement in nodes not fixed. Also, result is vector with same size
    // of reduced input
    virtual std::vector<Real> evaluate(std::vector<Real> uReduced) const override {
        // Pass reduced u to complete u, using fixedVarValues
        std::vector<Real> u = m_systemTransformations.reducedDisplacementToOriginalVector(uReduced);

        // Compute full result
        std::vector<Real> originalResult = m_function->evaluate(u);

        // Pass complete result to reduced result using reducedVarForVar
        std::vector<Real> reducedResult = m_systemTransformations.originalToReducedVector(originalResult);

        return reducedResult;
    }

    // Reduced form: meaning that u corresponds to displacement in nodes not fixed. Also, result is jacobian matrix
    // NOT for all nodes of the mesh
    virtual TripletMatrix<Triplet<Real>> jacobian(std::vector<Real> uReduced) const override {
        // Pass reduced u to complete u, using fixedVarValues
        std::vector<Real> u = m_systemTransformations.reducedDisplacementToOriginalVector(uReduced);

        // Compute full matrix result
        TripletMatrix<Triplet<Real>> originalResult = m_function->jacobian(u);

        // Pass complete result to reduced result using reducedVarForVar, eliminating rows and columns to be reduced
        TripletMatrix<Triplet<Real>> reducedResult = m_systemTransformations.originalToReducedMatrix(originalResult);

        return reducedResult;
    };

    virtual size_t size() const override {
        return m_size;
    }


private:
    std::shared_ptr<NonLinearElasticityFunction<Real>> m_function;
    SystemTransformations<Real> m_systemTransformations;
    size_t m_size;
};

#endif //NONLINEARELASTICITYFUNCTION_H

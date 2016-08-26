////////////////////////////////////////////////////////////////////////////////
// BaseCellOperations.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Organizes the computations we need to perform homogenization and shape
//      optimization on various types of microstructure base cells (triply
//      periodic, orthotropic, ...)
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  08/26/2016 03:08:14
////////////////////////////////////////////////////////////////////////////////
#ifndef BASECELLOPERATIONS_HH
#define BASECELLOPERATIONS_HH

#include <vector>
#include <memory>
#include <stdexcept>

#include <PeriodicHomogenization.hh>
#include <OrthotropicHomogenization.hh>
#include "SDConversions.hh"

enum class BaseCellType { TriplyPeriodic, Orthotropic };

template<class _Sim>
class BaseCellOperations {
public:
    static constexpr size_t N = _Sim::N;
    using VField  = typename _Sim::VField;
    using ETensor = typename _Sim::ETensor;

    static std::unique_ptr<BaseCellOperations<_Sim>> construct(BaseCellType type, _Sim &sim) {
        // It seems the actual construction must happen in a helper function
        // (_constructBaseCellOps)--couldn't declare this static member function
        // as a friend in derived classes possibly due to _Sim being a dependent
        // type...
        auto bco = _constructBaseCellOps(type, sim);
        bco->m_cellType = type;
        bco->m_solveCellProblems(sim);
        return bco;
    }

    virtual std::vector<typename _Sim::VField>
    solveAdjointCellProblems(std::vector<typename _Sim::VField> &adjointRHS) {
        if (adjointRHS.size() != flatLen(N)) throw std::runtime_error("Incorrect number of right-hand sides");

        std::vector<typename _Sim::VField> lambda;
        lambda.reserve(flatLen(N));
        for (size_t ij = 0; ij < flatLen(N); ++ij)
            lambda.emplace_back(m_solveProbeSystem(ij, adjointRHS[ij]));
        return lambda;
    }

    virtual ETensor homogenizedElaticityTensor() const = 0;

    const _Sim &sim() const { return m_sim; }
    const VField &w_ij(size_t ij)                         const { return m_w_ij.at(ij); }
    const std::vector<VField> &fluctuationDisplacements() const { return m_w_ij; }

    ////////////////////////////////////////////////////////////////////////////
    // Shape Derivative Conversions
    ////////////////////////////////////////////////////////////////////////////
    // TODO: make these non-virtual and manually dispatch based on BaseCellType
    // virtual diff_bdry_from_nsv_functional() const = 0;
    // virtual descent_from_diff_vol()         const = 0;
    // virtual diff_bdry_from_diff_vol()       const = 0;
    // virtual diff_bdry_from_nsv_functional() const = 0;

    virtual ~BaseCellOperations() { }

protected:
    BaseCellOperations(const _Sim &sim) : m_sim(sim) { }

    // Needs non-const sim
    virtual void m_solveCellProblems(_Sim &sim) = 0;
    virtual VField m_solveProbeSystem(size_t ij, const VField &rhs) const = 0;

    // Solved for by derived classes
    std::vector<typename _Sim::VField> m_w_ij;
    const _Sim &m_sim;

    BaseCellType m_cellType;
};

// Forward declare _constructBaseCellOps helper function so it can be friended.
template<class _Sim>
std::unique_ptr<BaseCellOperations<_Sim>> _constructBaseCellOps(BaseCellType type, _Sim &sim);

////////////////////////////////////////////////////////////////////////////////
// Triply Periodic Base Cell
////////////////////////////////////////////////////////////////////////////////
// For the traditional triply periodic homogenization, only a single system
// needs to be built, which is just stored in the simulator for historic
// reasons.
template<class _Sim>
class TriplyPeriodicBaseCellOperations : public BaseCellOperations<_Sim> {
public:
    static constexpr size_t N = _Sim::N;
    using Base    = BaseCellOperations<_Sim>;
    using VField  = typename Base::VField;
    using ETensor = typename Base::ETensor;

    virtual ETensor homogenizedElaticityTensor() const override {
        return PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(this->m_w_ij, this->m_sim);
    }

    virtual ~TriplyPeriodicBaseCellOperations() { }

protected:
    TriplyPeriodicBaseCellOperations(const _Sim &sim) : Base(sim) { }

    // Needs non-const sim
    virtual void m_solveCellProblems(_Sim &sim) override {
        PeriodicHomogenization::solveCellProblems(this->m_w_ij, sim, 1e-11);
    }

    virtual VField m_solveProbeSystem(size_t ij, const VField &rhs) const override {
        return this->m_sim.solve(rhs);
    }

    template<class _S2>
    friend std::unique_ptr<BaseCellOperations<_S2>> _constructBaseCellOps(BaseCellType type, _S2 &sim);
};

////////////////////////////////////////////////////////////////////////////////
// Orthotropic Base Cell
////////////////////////////////////////////////////////////////////////////////
// The different probing strains require different boundary conditions on the
// orthotropic base cell, meaning we need to store several probe systems.
//
// needs to be built, which is just stored in the simulator for historic
// reasons.
template<class _Sim>
class OrthotropicBaseCellOperations : public BaseCellOperations<_Sim> {
public:
    static constexpr size_t N = _Sim::N;
    using Base    = BaseCellOperations<_Sim>;
    using VField  = typename Base::VField;
    using ETensor = typename Base::ETensor;

    virtual ETensor homogenizedElaticityTensor() const override {
        return PeriodicHomogenization::Orthotropic::homogenizedElasticityTensorDisplacementForm(this->m_w_ij, this->m_sim);
    }

    virtual ~OrthotropicBaseCellOperations() { }

protected:
    OrthotropicBaseCellOperations(const _Sim &sim) : Base(sim) { }

    // Needs non-const sim
    virtual void m_solveCellProblems(_Sim &sim) override {
        m_probeSystems = PeriodicHomogenization::Orthotropic::solveCellProblems(this->m_w_ij, sim, 1e-11);
    }

    virtual VField m_solveProbeSystem(size_t ij, const VField &rhs) const override {
        return m_systemForProbe(ij).solve(rhs);
    }

    SPSDSystem<Real> &m_systemForProbe(size_t ij) const {
        if (ij >= flatLen(N)) throw std::runtime_error("probe condition index out of bounds");
        if (ij < N) return *m_probeSystems.at(0);
        else        return *m_probeSystems.at((ij - N) + 1);
    }

    std::vector<std::unique_ptr<SPSDSystem<Real>>> m_probeSystems;

    template<class _S2>
    friend std::unique_ptr<BaseCellOperations<_S2>> _constructBaseCellOps(BaseCellType type, _S2 &sim);
};

////////////////////////////////////////////////////////////////////////////////
// Helper method for allocating (but not initializing!) a base cell.
////////////////////////////////////////////////////////////////////////////////
template<class _Sim>
std::unique_ptr<BaseCellOperations<_Sim>>
_constructBaseCellOps(BaseCellType type, _Sim &sim) {
    std::unique_ptr<BaseCellOperations<_Sim>> bco;

    // Can't use make_unique since it would need to be a friend...
    if       (type == BaseCellType::TriplyPeriodic) bco = std::unique_ptr<TriplyPeriodicBaseCellOperations<_Sim>>(new TriplyPeriodicBaseCellOperations<_Sim>(sim));
    else if  (type == BaseCellType::Orthotropic   ) bco = std::unique_ptr<   OrthotropicBaseCellOperations<_Sim>>(new    OrthotropicBaseCellOperations<_Sim>(sim));
    else throw std::runtime_error("Unknown base cell type.");

    return bco;
}

#endif /* end of include guard: BASECELLOPERATIONS_HH */

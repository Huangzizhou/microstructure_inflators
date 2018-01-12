////////////////////////////////////////////////////////////////////////////////
// NonPeriodicCellOperations.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Organizes the computations needed to run shape optimization
//      Based on BaseCellOperations
*/
//  Author:  Davi Colli Tozoni (dctozoni), davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  1/9/18
////////////////////////////////////////////////////////////////////////////////

#ifndef NONPERIODICCELLOPERATIONS_H
#define NONPERIODICCELLOPERATIONS_H

#include <vector>
#include <memory>
#include <stdexcept>
#include <Future.hh>

#include <LinearElasticity.hh>
#include "SDConversions.hh"

using namespace std;

template<class _Sim>
class NonPeriodicCellOperations {
public:
    static constexpr size_t N = _Sim::N;
    using VField  = typename _Sim::VField;
    using ETensor = typename _Sim::ETensor;

    static std::unique_ptr<NonPeriodicCellOperations<_Sim>> construct(_Sim &sim, vector<CondPtr<N> > &bconds) {
        // It seems the actual construction must happen in a helper function
        // (_constructNonPeriodicCellOps)--couldn't declare this static member function
        // as a friend in derived classes possibly due to _Sim being a dependent
        // type...
        auto nco = _constructNonPeriodicCellOps(sim, bconds);
        nco->m_solveCellProblems(sim, bconds);
        return nco;
    }

    // Common laplace PDE with user defined boundary conditions
    virtual void m_solveCellProblems(_Sim &sim, vector<CondPtr<N> > &bconds) {
        VField f_zero(sim.numDoFs());
        f_zero.clear(); // to guarantee it is all zero

        sim.applyBoundaryConditions(bconds);
        sim.solve(f_zero);
    }

    // TODO: check what is the strong PDE form of the adjoint problem
    virtual VField m_solveProbeSystem(const VField &rhs) const {
        return this->m_sim.solve(rhs);
    }

    // Solve adjoint problem for the cell problem
    VField solveAdjointCellProblem(const VField &adjointRHS) const {
        return m_solveProbeSystem(adjointRHS);
    }

    const                _Sim & sim() const { return m_sim; }
    const typename _Sim::Mesh &mesh() const { return m_sim.mesh(); }
    const VField &displacement() const { return m_u; }

    ////////////////////////////////////////////////////////////////////////////
    // Shape Derivative Conversions
    // These are template functions and cannot be made virtual; we must handle
    // dispatch manually...
    ////////////////////////////////////////////////////////////////////////////
    template<size_t _SDDeg>
    ScalarOneForm<N> diff_bdry_from_nsv_functional(const std::vector<Interpolant<Real, N - 1, _SDDeg>> &sd) const {
        return SDConversions::diff_bdry_from_nsv_functional(sd, mesh());
    }

    template<typename T = Real>
    OneForm<T, N> diff_bdry_from_diff_vol(const OneForm<T, N> &diff_vol) const {
        return SDConversions::diff_bdry_from_diff_vol(diff_vol, sim());
    }

    VField descent_from_diff_bdry(const ScalarOneForm<N> &diff_bdry) const {
        return SDConversions::descent_from_diff_bdry(diff_bdry, sim());
    }

    VField descent_from_diff_vol(const ScalarOneForm<N> &diff_vol) const {
        return SDConversions::descent_from_diff_vol(diff_vol, sim());
    }

    virtual ~NonPeriodicCellOperations() { }

    NonPeriodicCellOperations(const _Sim &sim, vector<CondPtr<N> > &bconds) : m_sim(sim) { m_bconds = bconds; }

protected:

    const VField m_u;
    const _Sim &m_sim;
    vector<CondPtr<N> > m_bconds;
};

////////////////////////////////////////////////////////////////////////////////
// Helper method for allocating (but not initializing!) a cell.
////////////////////////////////////////////////////////////////////////////////
template<class _Sim>
std::unique_ptr<NonPeriodicCellOperations<_Sim>> _constructNonPeriodicCellOps(_Sim &sim, vector<CondPtr<_Sim::N> > &bconds) {
    std::unique_ptr<NonPeriodicCellOperations<_Sim>> nco;

    nco = std::unique_ptr<NonPeriodicCellOperations<_Sim>>(new NonPeriodicCellOperations<_Sim>(sim, bconds));

    return nco;
}

#endif //NONPERIODICCELLOPERATIONS_H


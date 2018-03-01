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

        //sim.applyNoRigidMotionConstraint();
        sim.applyBoundaryConditions(bconds);
        //m_u = sim.solve(f_zero);
        m_u = sim.solve();
    }

    // Solve adjoint problem for the cell problem
    VField solveAdjointCellProblem(const VField &adjointRHS) const {
        VField result;

        // Set all dirichlet boundary conditions to 0 (zero).
        m_sim.removeDirichletConditions();
        m_sim.removeNeumanConditions(); // TODO: necessary?!
        vector<CondPtr<N> > new_conds;
        for (unsigned i = 0; i < m_bconds.size(); i++) {
            CondPtr<N> cond = m_bconds[i];
            BoundaryCondition<N> *new_cond;
            if (const DirichletCondition<N> * dc = dynamic_cast<const DirichletCondition<N> *>(cond.get())) {
                VectorND<N> zero_vector;
                zero_vector.setZero();
                new_cond = new DirichletCondition<N>(dc->region, zero_vector, dc->componentMask);

                new_conds.push_back(CondPtr<N>(new_cond));
            }
        }

        // TODO verify: Neumann conditions are already embedded in adjointRHS!

        m_sim.applyBoundaryConditions(new_conds);
        result = m_sim.solve(adjointRHS);

        // Reapply old conditions
        m_sim.applyBoundaryConditions(m_bconds);

        return result;
    }

    // Change in the displacements due to mesh vertex perturbations delta_p
    VField deltaDisplacements(const VField &u, const VField &delta_p) const {
        VField result;

        // Set all dirichlet boundary conditions to 0 (zero).
        m_sim.removeDirichletConditions();
        vector<CondPtr<N> > new_conds;
        for (unsigned i = 0; i < m_bconds.size(); i++) {
            CondPtr<N> cond = m_bconds[i];
            BoundaryCondition<N> *new_cond;
            if (const DirichletCondition<N> * dc = dynamic_cast<const DirichletCondition<N> *>(cond.get())) {
                VectorND<N> zero_vector;
                zero_vector.setZero();
                new_cond = new DirichletCondition<N>(dc->region, zero_vector, dc->componentMask);

                new_conds.push_back(CondPtr<N>(new_cond));
            }
        }
        m_sim.applyBoundaryConditions(new_conds);

        auto rhs = -m_sim.applyDeltaStiffnessMatrix(u, delta_p);

        // Add part related to Neumann Boundary conditions
        VField deltaNeummanField = m_sim.deltaNeumannLoad(delta_p);
        rhs += deltaNeummanField;

        MSHFieldWriter writer("deltaDisplacements_validation", m_sim.mesh(), false);
        writer.addField("delta neumann field", deltaNeummanField);

        result = m_sim.solve(rhs);

        // Reapply old conditions
        m_sim.applyBoundaryConditions(m_bconds);

        return result;
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

    NonPeriodicCellOperations(_Sim &sim, vector<CondPtr<N> > &bconds) : m_sim(sim) { m_bconds = bconds; }

protected:

    VField m_u;
    _Sim &m_sim;
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


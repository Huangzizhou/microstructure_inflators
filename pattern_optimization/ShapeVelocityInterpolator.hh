////////////////////////////////////////////////////////////////////////////////
// ShapeVelocityInterpolator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      While the change in shape is determined entirely by the boundary
//      (normal) velocity field, it is often more accurate for the purposes of
//      discrete shape optimization to consider smooth perturbations of the
//      whole mesh instead of just the boundary. To do this, we interpolate the
//      boundary velocity by solving a Laplace equation.
//
//      To maintain periodicity, this interpolated shape velocity must be
//      periodic (but we don't need its boundary to stay clamped to the square
//      periodic cell). We achieve this by enforcing periodic boundary
//      conditions in the Laplace solve on the period cell.
//
//      Shape velocities of the "internal periodic vertices" (mesh boundary
//      vertices that would become internal after the mesh is periodically
//      tiled) are not read/applied as Dirichlet constraints (they should be
//      zero anyway, but the inflator is not required to guarantee this).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  03/19/2016 21:13:18
////////////////////////////////////////////////////////////////////////////////
#ifndef SHAPEVELOCITYINTERPOLATOR_HH
#define SHAPEVELOCITYINTERPOLATOR_HH

#include <Laplacian.hh>

class ShapeVelocityInterpolator {
public:
    template<class Sim>
    ShapeVelocityInterpolator(const Sim &sim) {
        const auto &mesh = sim.mesh();
        L = Laplacian::construct<1>(mesh);

        // Apply same periodic boundary conditions as sim.
        // Unfortunately, sim's periodic boundary conditions are on nodes, not
        // just vertices, so we have to reindex to get contiguous variables
        size_t nv = mesh.numVertices();
        m_varForVertex.reserve(nv);
        constexpr size_t NO_VAR = std::numeric_limits<size_t>::max();
        std::vector<size_t> varForDoF(sim.numDoFs(), NO_VAR);
        size_t numVars = 0;
        for (auto v : mesh.vertices()) {
            size_t d = sim.DoF(v.node().index());
            size_t &var = varForDoF.at(d);
            if (var == NO_VAR)
                var = numVars++;
            m_varForVertex.push_back(var);
        }

        // Apply periodic boundary conditions to L
        L.reindexVariables(numVars, m_varForVertex);

        // Determine which variables correspond to true boundary vertices (as
        // opposed to internal periodic vertices).
        // Vertices are true boundary vertices if any non-periodic boundary
        // element touches them.
        std::vector<bool> isTrueBoundaryVertex(mesh.numBoundaryVertices(), false);
        for (auto be : mesh.boundaryElements()) {
            if (be->isPeriodic) continue;
            for (size_t i = 0; i < be.numVertices(); ++i)
                isTrueBoundaryVertex.at(be.vertex(i).index()) = true;
        }

        // Determine vars corresponding to true boundary vertices (to be constrained)
        // Also, pick an arbitrary representative boundary vertex for a var
        // (for extracting Dirichlet condition)
        std::vector<bool> varFixed(numVars, false); // Fix each variable only once
        for (auto bv : mesh.boundaryVertices()) {
            if (isTrueBoundaryVertex.at(bv.index())) {
                size_t vari = m_varForVertex.at(bv.volumeVertex().index());
                if (varFixed.at(vari)) continue;
                m_bdryVars.push_back(vari);
                m_bdryVtxForBdryVar.push_back(bv.index());
                varFixed.at(vari) = true;
            }
        }
    }

    // Periodically interpolate per-boundary-vertex shape velocity vector field,
    // bdry_svel.
    // WARNING: non-periodic input fields will be made periodic by choosing the
    // value from an arbitrary one of the periodically identified vertices.
    template<class Sim>
    typename Sim::VField interpolate(const Sim &sim,
                                     const typename Sim::VField &bdry_svel) const {
        constexpr size_t N = Sim::N;
        const auto &mesh = sim.mesh();

        typename Sim::VField result;
        result.resizeDomain(mesh.numVertices()); // zero-inits

        std::vector<Real> bdryVarValues, zero(L.m), x;
        bdryVarValues.reserve(m_bdryVars.size());
        // Interpolate each component of the boundary velocity.
        for (size_t c = 0; c < N; ++c) {
            bdryVarValues.clear();
            // Extract boundary variable values from the per-
            for (size_t bvi : m_bdryVtxForBdryVar)
                bdryVarValues.push_back(bdry_svel(bvi)[c]);

            // TODO: allow changing the fixed variable constraint RHS without
            // rebuilding the system
            SPSDSystem<Real> Lsys(L);
            Lsys.fixVariables(m_bdryVars, bdryVarValues);
            Lsys.solve(zero, x);
            for (size_t vi = 0; vi < mesh.numVertices(); ++vi)
                result(vi)[c] = x.at(m_varForVertex.at(vi));
        }

        return result;
    }

private:
    // Periodic Laplacian
    TripletMatrix<> L;
    std::vector<size_t> m_varForVertex, m_bdryVars, m_bdryVtxForBdryVar;
};

#endif /* end of include guard: SHAPEVELOCITYINTERPOLATOR_HH */

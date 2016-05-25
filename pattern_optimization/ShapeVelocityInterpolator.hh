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

    // Find the discrete boundary velocity one-form, dvb, that satisfies
    //    dvb[vb] := dv[interpolate(vb)]
    // given the volume velocity one-form dv, for all periodic boundary
    // velocity fields vb.
    // In other words, applies [-L_bi L_ii^{-1} I] (the transpose of the matrix
    // applied for interpolation) to discrete volume velocity one-form dv.
    // Periodicity details:
    //    dv is non-periodic; it contains partial contributions on each
    //    identified periodic vertex (i.e. it gives the correct answer when
    //    dotted with a *periodic* nodal vector field).
    //    We introduce the following matrices mapping between vertex/dof field types:
    //       S^T: Sum identified vertices onto DoFs,                 S: copy DoFs to identified vertices
    //         B: Extract boundary DoFs from all DoFs,             B^T: distribute boundary DoFs to full DoF vector (zeros on internal DoFs)
    //         C: Extract internal DoFs from all DoFs,             C^T: distribute internal DoFs to full DoF vector (zeros on boundary DoFs)
    //         A: Average identified bdry vertices onto bdry DoFs, A^T: fractional distribution of bdry DoFs to identified bdry vertices
    //    dvb[vb] := dv[interpolate(vb)] = dv[S [-L_bi L_ii^-1 I]^T A vb]
    //             = lambda . vb
    //    ==> dvb = lambda = A^T [-L_bi L_ii^-1 I] S^T dv = A^T B  S^T dv - A^T L_bi L_ii^{-1} C S^T dv
    //                                                    = A^T B (S^T dv - L C^T    L_ii^{-1} C S^T dv)
    // WARNING, matrix A doesn't really "average" values from identified vertices:
    // For now, we just extract value from a single one of the identified vertices.
    // This is valid because we only operate on periodic fields, but it means
    // the resulting dvb will put zeros on all but one vertex in each identified vertex set.
    template<class Sim>
    typename Sim::VField adjoint(const Sim &sim,
                                 const typename Sim::VField &dv) const {
        assert(dv.domainSize() == sim.mesh().numVertices());

        // If there are no true boundary vertices, there is no adjoint velocity
        if (m_bdryVars.size() == 0) {
            typename Sim::VField result(sim.mesh().numBoundaryVertices());
            result.clear();
            return result;
        }
        SPSDSystem<Real> Lsys(L);

        // Fix boundary vars to zero
        // Now Lsys.solve() actually applies C^T L_ii^{-1} C
        Lsys.fixVariables(m_bdryVars, std::vector<Real>(m_bdryVars.size()));

        std::vector<Real> S_t_dv, C_t_Lii_inv_C_S_t_dv;
        typename Sim::VField dvb(sim.mesh().numBoundaryVertices());
        dvb.clear();
        for (size_t c = 0; c < Sim::N; ++c) {
            // Sum identified values onto the DoFs
            S_t_dv.assign(L.m, 0.0);
            for (size_t i = 0; i < m_varForVertex.size(); ++i)
                S_t_dv.at(m_varForVertex[i]) += dv(i)[c];

            Lsys.solve(S_t_dv, C_t_Lii_inv_C_S_t_dv); 
            auto dofValues = L.apply(C_t_Lii_inv_C_S_t_dv);
            // dofValues = S^T dv - L C^T    L_ii^{-1} C S^T dv
            for (size_t i = 0; i < dofValues.size(); ++i)
                dofValues[i] = S_t_dv[i] - dofValues[i];

            // Apply A^T B: extract boundary DoFs and then fractionally
            // distribute to bdry vertices (places full DoF value on a
            // single (true bdry) vtx per dof for now).
            for (size_t i = 0; i < m_bdryVars.size(); ++i)
                dvb(m_bdryVtxForBdryVar[i])[c] = dofValues.at(m_bdryVars[i]);
        }

        return dvb;
    }

private:
    // Periodic Laplacian
    TripletMatrix<> L;
    std::vector<size_t> m_varForVertex, m_bdryVars, m_bdryVtxForBdryVar;
};

#endif /* end of include guard: SHAPEVELOCITYINTERPOLATOR_HH */

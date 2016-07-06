#include "BoundaryPerturbationInflator.hh"

#include <UniformLaplacian.hh>
#include <MeshIO.hh>
#include <Future.hh>

template<size_t N>
BoundaryPerturbationInflator<N>::BoundaryPerturbationInflator(
        const std::string &meshPath, Real epsilon)
{
    std::vector<MeshIO::IOVertex>  inVertices;
    std::vector<MeshIO::IOElement> inElements;
    auto type = MeshIO::load(meshPath, inVertices, inElements);

    size_t dim;
    if      (type == MeshIO::MESH_TET) dim = 3;
    else if (type == MeshIO::MESH_TRI) dim = 2;
    else    throw std::runtime_error("Mesh must be triangle or tet.");
    
    if (dim != N) throw std::runtime_error("Mesh/inflator dimension match.");

    m_setMesh(inVertices, inElements, epsilon);
}

template<size_t N>
BoundaryPerturbationInflator<N>::BoundaryPerturbationInflator(
        const std::vector<MeshIO::IOVertex>  &inVertices,
        const std::vector<MeshIO::IOElement> &inElements, Real epsilon)
{
    m_setMesh(inVertices, inElements, epsilon);
}

////////////////////////////////////////////////////////////////////////////
// Inflation
////////////////////////////////////////////////////////////////////////////
// Solve for all coordinate variables using uniform Laplacian, with params
// and bounding box positions as Dirichlet constraints
// Note, "params" are treated as *offsets* to m_origParams
template<size_t N>
void BoundaryPerturbationInflator<N>::inflate(const std::vector<Real> &params)
{
    assert(params.size() == m_numParams);

    std::vector<size_t> fixedVars;
    std::vector<Real> fixedVarValues;
    // conservative allocation: m_numParams/N params per dim on average
    fixedVars.reserve(2 * N + m_numParams);
    fixedVarValues.reserve(2 * N + m_numParams);

    const auto bbox = m_mesh->boundingBox();

    m_vertices.clear(), m_vertices.resize(m_mesh->numVertices());
    m_elements.clear(), m_elements.assign(m_mesh->numElements(), MeshIO::IOElement(N + 1));
    for (auto e : m_mesh->elements()) {
        for (size_t c = 0; c < N + 1; ++c)
            m_elements[e.index()][c] = e.vertex(c).index();
    }

    if (m_noPerturb) {
        for (Real p : params) {
            if (p != 0.0)
                throw std::runtime_error("Requested no-perturb mode with nonzero parameter vector");
        }
        for (auto v : m_mesh->vertices())
            m_vertices[v.index()] = v.node()->p;
        return;
    }

    for (size_t d = 0; d < N; ++d) {
        SPSDSystem<Real> L;
        UniformLaplacian::assemble(*m_mesh, L, m_varForCoordinate[d]);

        fixedVars.clear(); fixedVarValues.clear();
        // Add the special var Dirichlet constraints
        for (size_t i = 0; i < 2 * N; ++i) {
            fixedVars.push_back(i);
            fixedVarValues.push_back(i < N ? bbox.minCorner[i]
                                           : bbox.maxCorner[i - N]);
        }
        // Add the parameter Dirichlet constraints
        for (size_t vari = 0; vari < m_numVars[d]; ++vari) {
            size_t p = m_paramForVariable[d][vari];
            if (p != NONE) {
                fixedVars.push_back(vari);
                fixedVarValues.push_back(params.at(p) + m_origParams[p]);
            }
        }
        
        // Solve for all variables
        L.fixVariables(fixedVars, fixedVarValues);
        std::vector<Real> zero(m_numVars[d], 0.0), x;
        L.solve(zero, x);
        assert(x.size() == m_numVars[d]);

        // Read off vertex coordinates
        for (size_t vi = 0; vi < m_mesh->numVertices(); ++vi)
            m_vertices.at(vi)[d] = x.at(m_varForCoordinate[d].at(vi));
    }

    // Reposition vertices into our mesh data structure to support normal
    // and shape velocity computations
    m_mesh->setNodePositions(m_vertices);
}


// Read off the parameter values from a particular per-boundary-vertex
// vector field, verifying its consistency with the periodic constraints
// If guaranteeing consistent boundary vertex enumerations across multiple
// FEMMesh instances becomes a problem, we could change this to take a
// per-volume-vertex field.
template<size_t N>
ScalarField<Real> BoundaryPerturbationInflator<N>::
paramsFromBoundaryVField(const VectorField<Real, N> &values) const
{
    ScalarField<Real> result(m_numParams);
    std::vector<bool> isSet(m_numParams, false);
    assert(values.domainSize() == m_mesh->numBoundaryVertices());
    for (auto bv : m_mesh->boundaryVertices()) {
        for (size_t d = 0; d < N; ++d) {
            Real val = values(bv.index())[d];
            size_t var = m_varForCoordinate[d].at(bv.volumeVertex().index());
            size_t p = m_paramForVariable[d].at(var);
            if (p != NONE) {
                if (!isSet.at(p)) {
                    result[p] = val;
                    isSet.at(p) = true;
                }
                if (std::abs(result[p] - val) > 1e-6)
                    std::cerr << "WARNING: boundary values violate periodicity constraints" << std::endl;;
            }
        }
    }
    return result;
}

// Get the boundary vector field corresponding to "params" (i.e. the
// inverse of paramsFromBoundaryVField)
template<size_t N>
VectorField<Real, N> BoundaryPerturbationInflator<N>::
boundaryVFieldFromParams(const ScalarField<Real> &params) const
{
    VectorField<Real, N> result(m_mesh->numBoundaryVertices());
    result.clear();
    for (auto bv : m_mesh->boundaryVertices()) {
        for (size_t d = 0; d < N; ++d) {
            auto var = m_varForCoordinate[d].at(bv.volumeVertex().index());
            size_t p = m_paramForVariable[d].at(var);
            if (p != NONE) result(bv.index())[d] = params[p];
        }
    }

    return result;
}

// Initialize the boundary perturbation inflator for a particular mesh.
template<size_t N>
void BoundaryPerturbationInflator<N>::m_setMesh(
    const std::vector<MeshIO::IOVertex>  &inVertices,
    const std::vector<MeshIO::IOElement> &inElements, Real epsilon)
{
    m_mesh = Future::make_unique<Mesh>(inElements, inVertices);

    // Note: pc matches every node, not just vertex-collocated nodes!
    // But vertex nodes are a prefix of all nodes, so we can ignore this.
    PeriodicCondition<N> pc(*m_mesh, epsilon);
    ////////////////////////////////////////////////////////////////////////
    // Determine variables (apply periodic coordinate constraints)
    ////////////////////////////////////////////////////////////////////////
    m_numVars.fill(2 * N); // variables 0..2N-1 always store the periodic
                           // face coordinates
    m_varForCoordinate.fill(std::vector<size_t>(m_mesh->numVertices(),
                                                size_t(NONE)));

    // Func creating new component d var. NOTE & captures this, not members
    auto createVar = [&](size_t d) { return m_numVars[d]++; };

    // All interior coordinates get distinct variables.
    for (auto v : m_mesh->vertices()) {
        if (v.isBoundary()) continue;
        for (size_t d = 0; d < N; ++d)
            m_varForCoordinate[d].at(v.index()) = createVar(d);
    }
    std::vector<size_t> presentPairs;
    for (auto bv : m_mesh->boundaryVertices()) {
        size_t bvi = bv.index();
        size_t vvi = bv.volumeVertex().index();

        for (size_t d = 0; d < N; ++d) {
            // Link vertices at the min/max faces to the special variables
            // holding min/max bbox coordinates
            int minMax = pc.bdryNodeOnMinOrMaxPeriodCellFace(bvi, d);
            if      (minMax == -1) m_varForCoordinate[d].at(vvi) = d;
            else if (minMax ==  1) m_varForCoordinate[d].at(vvi) = N + d;
            else {
                // For the coordinates not clamped to min/max, introduce new
                // variables shared by all identified nodes.
                const size_t currVar = m_varForCoordinate[d].at(vvi);
                const size_t newVar = (currVar != NONE) ? currVar
                                                        : createVar(d);
                for (size_t iv : pc.identifiedNodes(vvi)) {
                    assert(pc.bdryNodeOnMinOrMaxPeriodCellFace(m_mesh->vertex(iv).boundaryVertex().index(), d) == 0);
                    assert(m_varForCoordinate[d].at(iv) == currVar);
                    m_varForCoordinate[d].at(iv) = newVar;
                }
            }
        }
    }

    // Sanity check--all vertices should be linked to a variable.
    for (size_t vvi = 0; vvi < m_mesh->numVertices(); ++vvi) {
        for (size_t d = 0; d < N; ++d)
            assert(m_varForCoordinate[d].at(vvi) != NONE);
    }

    ////////////////////////////////////////////////////////////////////////
    // Determine parameters (find true microstructure boundary vertices)
    ////////////////////////////////////////////////////////////////////////
    // Initialize with no parameters.
    for (size_t d = 0; d < N; ++d)
        m_paramForVariable[d].assign(m_numVars[d], size_t(NONE));
    m_numParams = 0;
    // Func creating a new parameter. NOTE & captures this, not members
    auto createParam = [&]() { return m_numParams++; };
    // True boundary vertices are those with at least one incident
    // non-periodic-boundary element.
    std::vector<bool> isTrueBdryVertex(m_mesh->numBoundaryVertices(), false);
    for (auto be : m_mesh->boundaryElements()) {
        if (pc.isPeriodicBE(be.index())) continue;
        for (size_t c = 0; c < be.numVertices(); ++c)
            isTrueBdryVertex.at(be.vertex(c).index()) = true;
    }
    // Create a new param for every variable of each true boundary vertex
    // (But do not create parameters for "special variables" 0..2N-1 that
    //  are constrained to be the period cell min/max coordinates.)
    for (auto bv : m_mesh->boundaryVertices()) {
        if (!isTrueBdryVertex.at(bv.index())) continue;
        size_t vvi = bv.volumeVertex().index();
        size_t numAssigned = 0;
        for (size_t d = 0; d < N; ++d) {
            size_t var = m_varForCoordinate[d].at(vvi);
            if (var >= 2 * N) {
                if (m_paramForVariable[d].at(var) == NONE)
                    m_paramForVariable[d].at(var) = createParam();
                ++numAssigned;
            }
        }
        // True boundary vertices should have at least one associated
        // variable (assuming tiled pattern is manifold)
        assert(numAssigned > 0); 
    }

    // Get the original, pre-perturbation parameter values
    VectorField<Real, N> bdryPositions(m_mesh->numBoundaryVertices());
    for (auto bv : m_mesh->boundaryVertices())
        bdryPositions(bv.index()) = bv.volumeVertex().node()->p;
    m_origParams = paramsFromBoundaryVField(bdryPositions);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: 2D and 3D inflators.
////////////////////////////////////////////////////////////////////////////////
template class BoundaryPerturbationInflator<2>;
template class BoundaryPerturbationInflator<3>;
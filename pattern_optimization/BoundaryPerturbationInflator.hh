////////////////////////////////////////////////////////////////////////////////
// BoundaryPerturbationInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Inflator taking an initial mesh and using the boundary vertex position
//      offsets as parameters.
//
//      Vertices on the periodic boundary are alowed to move, but are
//      constrained to stay on their original periodic boundaries (i.e. cannot
//      leave or enter a new one--that would mean a topology change). Also,
//      they are constrained to move consistently with their periodic pair(s).
//      WARNING: vertices that start on a corner are constrained to stay there!
//
//      Only the positions of the "true" boundary vertices are specified as
//      parameters. All other position variables are solved for using a uniform
//      graph Laplacian. "True" boundary vertices are those with at least one
//      non-periodic boundary element incident. In other words, they are
//      the vertices that are still on the boundary after the period cell's
//      identified faces have been stitched together. Note that the true
//      boundary vertices lying on P periodic boundaries (0<=P<=N) will only
//      have N-P parameters due to the periodicity constraints, and these will
//      be shared by all 2^P identified vertices.
//      
//      First, periodicity constraints are enforced by constructing a reduced
//      set of variables. Then then we construct a set of parameters (one
//      parameter for each "true" boundary variable) from the variable set. The
//      i^th coordinate of each vertex is expressed in terms of variable
//      indices by array m_varForCoordinate[i], and the parameter corresponding
//      to each variable is given by m_paramForVariable[i] (which is NONE for
//      dependent variables).
//
//      Given a set of parameters, vertex coordinates are computed by:
//          For each coordinate i in 0 to N-1
//              1) Solve uniform Laplacian with parameter values (and periodic
//                 face coordinates) as Dirichlet constraints.
//              2) Decode the variables into vertex coordinates using
//                 m_paramForVariable[i]
//
//      Given a set of vertex coordinates/offsets, parameter values/offsets are
//      computed by:
//      For each coordinate i in 0 to N-1
//          1) Loop over vertices and extract variable values using
//             m_varForCoordinate[i]. Validate consistent extraction (i.e.
//             input coordinates should respectd periodic constraints)
//          2) Loop over variables and extract parameters using
//             m_paramForVariable[i]
//
//      There are 2 * N special variables, indexed 0..2N - 1 that hold the
//      periodic face coordinates:
//      0: min x, 1: min y[, 2: min z], N: max x, N + 1: max y[, N + 2: max z]
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/19/2015 02:14:38
////////////////////////////////////////////////////////////////////////////////
#ifndef BOUNDARYPERTURBATIONINFLATOR_HH
#define BOUNDARYPERTURBATIONINFLATOR_HH
#include <vector>
#include <array>

#include <GlobalBenchmark.hh>

#include "InflatorBase.hh"
#include <UniformLaplacian.hh>
#include <FEMMesh.hh>
#include <Fields.hh>
#include <BoundaryConditions.hh>
#include <GaussQuadrature.hh>

template<size_t N>
class BoundaryPerturbationInflator : public InflatorBase<BoundaryPerturbationInflator<N>> {
public:
    using Base = InflatorBase<BoundaryPerturbationInflator<N>>;
    using Base::m_vertices;
    using Base::m_elements;
    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    using Mesh = FEMMesh<N, 1, VectorND<N>>;
    using NSV = NormalShapeVelocity<N>;

    BoundaryPerturbationInflator(const std::vector<MeshIO::IOVertex>  &inVertices,
                                 const std::vector<MeshIO::IOElement> &inElements)
        : m_mesh(inElements, inVertices), m_pc(m_mesh)
    {
        ////////////////////////////////////////////////////////////////////////
        // Determine variables (apply periodic coordinate constraints)
        ////////////////////////////////////////////////////////////////////////
        m_numVars.fill(2 * N); // variables 0..2N-1 always store the periodic
                               // face coordinates
        m_varForCoordinate.fill(std::vector<size_t>(m_mesh.numVertices(),
                                                    size_t(NONE)));

        // Func creating new component d var. NOTE & captures this, not members
        auto createVar = [&](size_t d) { return m_numVars[d]++; };

        // All interior coordinates get distinct variables.
        for (auto v : m_mesh.vertices()) {
            if (v.isBoundary()) continue;
            for (size_t d = 0; d < N; ++d)
                m_varForCoordinate[d].at(v.index()) = createVar(d);
        }
        std::vector<size_t> presentPairs;
        for (auto bv : m_mesh.boundaryVertices()) {
            size_t bvi = bv.index();
            size_t vvi = bv.volumeVertex().index();

            for (size_t d = 0; d < N; ++d) {
                // Link vertices at the min/max faces to the special variables
                // holding min/max bbox coordinates
                int minMax = m_pc.bdryNodeOnMinOrMaxPeriodCellFace(bvi, d);
                if      (minMax == -1) m_varForCoordinate[d].at(vvi) = d;
                else if (minMax ==  1) m_varForCoordinate[d].at(vvi) = N + d;
                else {
                    // For the coordinates not clamped to min/max, introduce new
                    // variables shared by all identified nodes.
                    const size_t currVar = m_varForCoordinate[d].at(vvi);
                    const size_t newVar = (currVar != NONE) ? currVar
                                                            : createVar(d);
                    for (size_t iv : m_pc.identifiedNodes(vvi)) {
                        assert(m_pc.bdryNodeOnMinOrMaxPeriodCellFace(m_mesh.vertex(iv).boundaryVertex().index(), d) == 0);
                        assert(m_varForCoordinate[d].at(iv) == currVar);
                        m_varForCoordinate[d].at(iv) = newVar;
                    }
                }
            }
        }

        // Sanity check--all vertices should be linked to a variable.
        for (size_t vvi = 0; vvi < m_mesh.numVertices(); ++vvi) {
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
        std::vector<bool> isTrueBdryVertex(m_mesh.numBoundaryVertices(), false);
        for (auto be : m_mesh.boundaryElements()) {
            if (m_pc.isPeriodicBE(be.index())) continue;
            for (size_t c = 0; c < be.numVertices(); ++c)
                isTrueBdryVertex.at(be.vertex(c).index()) = true;
        }
        // Create a new param for every variable of each true boundary vertex
        // (But do not create parameters for "special variables" 0..2N-1 that
        //  are constrained to be the period cell min/max coordinates.)
        for (auto bv : m_mesh.boundaryVertices()) {
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
            // True boundary vertices should have at least one assocated
            // variable (assuming tiled pattern is manifold)
            assert(numAssigned > 0); 
        }

        // Get the original, pre-perturbation parameter values
        std::vector<VectorND<N>> bdryPositions;
        for (auto bv : m_mesh.boundaryVertices())
            bdryPositions.push_back(bv.volumeVertex().node()->p);
        m_origParams = extractParamsFromBoundaryValues(bdryPositions);
    }

    void setNoPerturb(bool noPerturb) { m_noPerturb = noPerturb; }
        
    // Solve for all coordinate variables using uniform Laplacian, with params
    // and bounding box positions as Dirichlet constraints
    // Note, "params" are treated as *offsets* to m_origParams
    void inflate(const std::vector<Real> &params) {
        assert(params.size() == m_numParams);

        std::vector<size_t> fixedVars;
        std::vector<Real> fixedVarValues;
        // conservative allocation: m_numParams/N params per dim on average
        fixedVars.reserve(2 * N + m_numParams);
        fixedVarValues.reserve(2 * N + m_numParams);

        const auto bbox = m_mesh.boundingBox();

        m_vertices.clear(), m_vertices.resize(m_mesh.numVertices());
        m_elements.clear(), m_elements.assign(m_mesh.numElements(), MeshIO::IOElement(N + 1));
        for (auto e : m_mesh.elements()) {
            for (size_t c = 0; c < N + 1; ++c)
                m_elements[e.index()][c] = e.vertex(c).index();
        }

        if (m_noPerturb) {
            for (Real p : params) {
                if (p != 0.0)
                    throw std::runtime_error("Requested no-perturb mode with nonzero parameter vector");
            }
            for (auto v : m_mesh.vertices())
                m_vertices[v.index()] = v.node()->p;
            return;
        }

        for (size_t d = 0; d < N; ++d) {
            SPSDSystem<Real> L;
            UniformLaplacian::assemble(m_mesh, L, m_varForCoordinate[d]);

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
                    fixedVarValues.push_back(params.at(p) + m_origParams.at(p));
                }
            }
            
            // Solve for all variables
            L.fixVariables(fixedVars, fixedVarValues);
            std::vector<Real> zero(m_numVars[d], 0.0), x;
            L.solve(zero, x);
            assert(x.size() == m_numVars[d]);

            // Read off vertex coordinates
            for (size_t vi = 0; vi < m_mesh.numVertices(); ++vi)
                m_vertices.at(vi)[d] = x.at(m_varForCoordinate[d].at(vi));
        }

        // Reposition vertices into our mesh data structure to support normal
        // and shape velocity computations
        m_mesh.setNodePositions(m_vertices);
    }

    // Read off the parameter values from a particular per-boundary-vertex
    // vector field, verifying its consistency with the periodic constraints
    // If guaranteeing consistent boundary vertex enumerations across multiple
    // FEMMesh instances becomes a problem, we could change this to take a
    // per-volume-vertex field.
    std::vector<Real> extractParamsFromBoundaryValues(
            const std::vector<VectorND<N>> &values) const {
        std::vector<Real> result(m_numParams);
        std::vector<bool> isSet(m_numParams, false);
        assert(values.size() == m_mesh.numBoundaryVertices());
        for (auto bv : m_mesh.boundaryVertices()) {
            for (size_t d = 0; d < N; ++d) {
                Real val = values[bv.index()][d];
                size_t var = m_varForCoordinate[d].at(bv.volumeVertex().index());
                size_t p = m_paramForVariable[d].at(var);
                if (p != NONE) {
                    if (!isSet.at(p)) {
                        result.at(p) = val;
                        isSet.at(p) = true;
                    }
                    if (std::abs(result.at(p) - val) > 1e-6)
                        std::cerr << "WARNING: boundary values violate periodicity constriants";
                }
            }
        }
        return result;
    }

    // Boundary vertex normal vector field (0 on periodic boundary) Uses area
    // for averaging.
    VectorField<Real, N> normals() const {
        VectorField<Real, N> result(m_mesh.numBoundaryVertices());
        result.clear();
        for (auto be : m_mesh.boundaryElements()) {
            if (m_pc.isPeriodicBE(be.index())) continue;
            for (size_t c = 0; c < be.numVertices(); ++c)
                result(be.vertex(c).index()) += be->volume() * be->normal();
        }
        for (size_t i = 0; i < result.domainSize(); ++i) {
            Real norm = result(i).norm();
            if (norm > 1e-6)
                result(i) /= norm;
        }
        return result;
    }

    ParameterType parameterType(size_t /* p */) const {
        return ParameterType::Offset;
    }

    bool isPrintable(const ScalarField<Real> &/* p */) const {
        // TODO
        return true;
    }

    template<class Mesh2>
    std::vector<NSV> computeShapeNormalVelocities(const Mesh2 &/* mesh */) const {
        throw std::runtime_error("Normal shape velocity unsupported for BoundaryPerturbationInflator");
    }

    // Takes a boundary interpolant representing the shape derivative of a
    // scalar objective function and computes the gradient with respect the
    // boundary coordinate parameters
    // WARNING: this is not actually a gradient, but rather a
    // steepest descent direction with respect to an approximate surface
    // distance metric.
    // This means it probably WON'T work with fancier optimization methods.
    template<class SDInterp>
    ScalarField<Real> gradientFromShapeDerivative(const std::vector<SDInterp> &sd) const {
        // BENCHMARK_START_TIMER("gradientFromShapeDerivative");
        assert(sd.size() == m_mesh.numBoundaryElements());

        ScalarField<Real> gradp(numParameters());
        gradp.clear();
        using NSVI = typename NSV::value_type;
        NSVI nsv;
        // Total area of all boundary elements adjacent to the vertex a
        // parameter controls (after conceptually stitching together the
        // identified periodic cell faces).
        std::vector<Real> paramArea(numParameters(), 0.0);
        for (auto be : m_mesh.boundaryElements()) {
            if (m_pc.isPeriodicBE(be.index())) continue; // periodic boundary elements don't count...
            auto normal = be->normal();
            const auto &sd_be = sd.at(be.index());
            for (size_t v = 0; v < be.numVertices(); ++v) {
                static_assert(be.numVertices() == nsv.size(), "Boundary element and NSV size mismatch");
                auto bv = be.vertex(v);
                size_t vvi = bv.volumeVertex().index();
                for (size_t d = 0; d < N; ++d) {
                    size_t var   = m_varForCoordinate[d].at(vvi);
                    size_t p     = m_paramForVariable[d].at(var);
                    if (p != NONE) {
                        nsv = 0;
                        nsv[v] = normal[d]; // v_p dot normal = e_d dot normal
                        gradp[p] += Quadrature<N - 1, NSVI::Deg + SDInterp::Deg>::
                            integrate([&] (const VectorND<be.numVertices()> &pt) {
                                return nsv(pt) * sd_be(pt);
                            }, be->volume());
                        paramArea[p] += be->volume();
                    }
                }
            }
        }

        // Normalize by parameter area (to make shape velocity less
        // mesh-dependent)
        for (size_t p = 0; p < numParameters(); ++p)
            gradp[p] /= paramArea[p];

        // BENCHMARK_STOP_TIMER("gradientFromShapeDerivative");
        return gradp;
    }

    // Takes a boundary interpolant representing the shape derivative of a
    // scalar objective function and computes the gradient with respect the
    // boundary coordinate parameters
    template<class SDInterp>
    ScalarField<Real> gradientFromShapeDerivativeVertexNormalVersion(const
            std::vector<SDInterp> &sd) const {
        assert(sd.size() == m_mesh.numBoundaryElements());

        VectorField<Real, N> vtxNormal = normals();
        const size_t nbv = m_mesh.numBoundaryVertices();
        assert(vtxNormal.domainSize() == nbv);

        // Unweighted average of steepest steepest ascent normal shape velocity
        // for each vertex, multiplied by "vertex area" (1/(N+1) of adjacent N-simplex area) 
        // (In the future, could integrate over voronoi region)
        ScalarField<Real> totalVertexNSV(nbv);
        totalVertexNSV.clear();
        std::vector<int> incidentCount(nbv, 0);
        ScalarField<Real> vertexArea(nbv);
        vertexArea.clear();
        for (auto be : m_mesh.boundaryElements()) {
            const auto &sd_be = sd.at(be.index());
            static_assert((SDInterp::Deg == 0) || (SDInterp::numNodalValues >= be.numVertices()), "Invalid SDInterp size");
            for (size_t v = 0; v < be.numVertices(); ++v) {
                size_t bvi = be.vertex(v).index();
                Real val = sd_be[0]; // deg 0 case
                if (SDInterp::Deg > 0)
                    val = sd_be[v];
                totalVertexNSV[bvi] += val;
                vertexArea[bvi] += be->volume() / be.numVertices();
                ++incidentCount.at(bvi);
            }
        }
        std::vector<VectorND<N>> perturbation(nbv);
        for (size_t bvi = 0; bvi < nbv; ++bvi) {
            Real nsv = totalVertexNSV[bvi] / Real(incidentCount.at(bvi));
            perturbation[bvi] = vtxNormal(bvi) * (nsv * vertexArea[bvi]);
        }

        return ScalarField<Real>(extractParamsFromBoundaryValues(perturbation));
    }

    // Get the boundary vector field corresponding to "params" (i.e. the
    // inverse of extractParamsFromBoundaryValues)
    VectorField<Real, N> boundaryVectorField(const ScalarField<Real> &params) const {
        VectorField<Real, N> result(m_mesh.numBoundaryVertices());
        result.clear();
        for (auto bv : m_mesh.boundaryVertices()) {
            for (size_t d = 0; d < N; ++d) {
                auto var = m_varForCoordinate[d].at(bv.volumeVertex().index());
                size_t p = m_paramForVariable[d].at(var);
                if (p != NONE) result(bv.index())[d] = params[p];
            }
        }

        return result;
    }

    const Mesh &mesh() const { return m_mesh; }
    const PeriodicCondition<N> &pc() const  { return m_pc; }

    size_t numParameters() const { return m_numParams; }

private:
    // Vector of indices into the coordinate variables
    std::array<std::vector<size_t>, N> m_varForCoordinate;
    std::array<std::vector<size_t>, N> m_paramForVariable;

    size_t m_numParams;
    std::array<size_t, N> m_numVars;

    // The original boundary coordinates pre-perturbation (i.e., the
    // coordinates to which offset parameters are applied)
    std::vector<Real> m_origParams;

    Mesh m_mesh;

    // Avoid perturbing interior vertices when parameter vector is zero.
    // (Useful for remeshing gradient descent.)
    bool m_noPerturb = false;

    // Note: pc matches every node, not just vertex nodes!
    // But vertex nodes are a prefix of all nodes, so we can ignore this.
    PeriodicCondition<N> m_pc;
};

#endif /* end of include guard: BOUNDARYPERTURBATIONINFLATOR_HH */

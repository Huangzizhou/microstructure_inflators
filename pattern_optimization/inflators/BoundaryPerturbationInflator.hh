////////////////////////////////////////////////////////////////////////////////
// BoundaryPerturbationInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Inflator taking an initial mesh and using the boundary vertex position
//      offsets as parameters.
//
//      Vertices on the periodic boundary are allowed to move, but are
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
//
//  Alternate versions:
//    - Leave internal vertices unperturbed.
//    - Solve for internal vertex perturbation using Laplacian w/ Dirichlet
//      boundary conditions.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/19/2015 02:14:38
////////////////////////////////////////////////////////////////////////////////
#ifndef BOUNDARYPERTURBATIONINFLATOR_HH
#define BOUNDARYPERTURBATIONINFLATOR_HH
#include <vector>
#include <array>
#include <limits>

#include <GlobalBenchmark.hh>

#include "../Inflator.hh"
#include <UniformLaplacian.hh>
#include <FEMMesh.hh>
#include <Fields.hh>
#include <BoundaryConditions.hh>

template<size_t N>
class BoundaryPerturbationInflator : public Inflator<N> {
public:
    using Base = Inflator<N>;
    using Base::m_vertices;
    using Base::m_elements;

    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    using Mesh = FEMMesh<N, 1, VectorND<N>>;

    BoundaryPerturbationInflator(const std::string &meshPath,
                                 Real epsilon = 1e-5);

    BoundaryPerturbationInflator(const std::vector<MeshIO::IOVertex>  &inVertices,
                                 const std::vector<MeshIO::IOElement> &inElements,
                                 Real epsilon = 1e-5);

    // Non periodic case
    BoundaryPerturbationInflator(const std::vector<MeshIO::IOVertex>  &inVertices,
                                 const std::vector<MeshIO::IOElement> &inElements,
                                 std::vector<CondPtr<N> > &bconds,
                                 Real epsilon = 1e-5);

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override;
public:

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity
    ////////////////////////////////////////////////////////////////////////////
    // Read off the parameter values from a particular per-boundary-vertex
    // vector field, verifying its consistency with the periodic constraints
    // If guaranteeing consistent boundary vertex enumerations across multiple
    // FEMMesh instances becomes a problem, we could change this to take a
    // per-volume-vertex field.
    virtual ScalarField<Real> paramsFromBoundaryVField(const VectorField<Real, N> &values) const override;

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return false; }
    virtual size_t numParameters() const override { return m_numParams; }
    virtual std::vector<Real> defaultParameters() const override { return std::vector<Real>(m_numParams); }
    virtual ParameterType parameterType(size_t /* p */) const override {
        return ParameterType::Offset;
    }
    virtual bool isPrintable(const std::vector<Real> &/* p */) override {
        // TODO
        return true;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Boundary perturbation-specific
    ////////////////////////////////////////////////////////////////////////////
    void setNoPerturb(bool noPerturb) { m_noPerturb = noPerturb; }

    // Boundary vertex normal vector field (0 on periodic boundary) Uses area
    // for averaging.
    VectorField<Real, N> boundaryVertexNormals() const {
        VectorField<Real, N> result(m_mesh->numBoundaryVertices());
        result.clear();
        for (auto be : m_mesh->boundaryElements()) {
            if (m_isPeriodicBE.at(be.index())) continue;
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

    // Get the boundary vector field corresponding to "params" (i.e. the
    // inverse of paramsFromBoundaryVField)
    virtual VectorField<Real, N> boundaryVFieldFromParams(const ScalarField<Real> &params) const override;

    const Mesh &mesh() const { return *m_mesh; }

    bool isInsideBoundaryCondition(size_t vi, std::vector<CondPtr<N> > &bconds) const {
        for (CondPtr<N> bc : bconds) {
            if (bc->containsPoint(m_mesh->vertex(vi).node()->p))
                return true;
        }

        return false;
    }

    virtual ~BoundaryPerturbationInflator() { }

private:
    // Vector of indices into the coordinate variables
    std::array<std::vector<size_t>, N> m_varForCoordinate;
    std::array<std::vector<size_t>, N> m_paramForVariable;
    std::vector<bool> m_isPeriodicBE;
    std::array<std::vector<bool>, N> m_bcVertexVariable; // say if variable is from boundary condition node
    std::array<std::vector<Real>, N> m_bcVertexValue; // saves variable constant value when

    size_t m_numParams;
    std::array<size_t, N> m_numVars;

    // The original boundary coordinates pre-perturbation (i.e., the
    // coordinates to which offset parameters are applied)
    ScalarField<Real> m_origParams;

    // Avoid perturbing interior vertices when parameter vector is zero.
    // (Useful for remeshing gradient descent.)
    bool m_noPerturb = false;
    bool m_isPeriodicMesh = true;

    std::unique_ptr<Mesh> m_mesh;
    std::vector<CondPtr<N> > m_bconds;

    void m_setMesh(const std::vector<MeshIO::IOVertex>  &inVertices,
                   const std::vector<MeshIO::IOElement> &inElements, Real epsilon);

    // Non periodic case
    void m_setMesh(const std::vector<MeshIO::IOVertex>  &inVertices,
                   const std::vector<MeshIO::IOElement> &inElements,
                   std::vector<CondPtr<N> > bconds,
                   Real epsilon);
};

#endif /* end of include guard: BOUNDARYPERTURBATIONINFLATOR_HH */

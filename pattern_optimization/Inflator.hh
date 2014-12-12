////////////////////////////////////////////////////////////////////////////////
// Inflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Wrapper providing a unified interface for the 2D and 3D wire inflators.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/11/2014 04:52:34
////////////////////////////////////////////////////////////////////////////////
#ifndef INFLATOR_HH
#define INFLATOR_HH
#include <WireInflator2D.h>
#include <Wires/Interfaces/PeriodicExploration.h>

enum class ParameterType { Thickness, Offset };

template<class _Derived>
class InflatorBase {
public:
    void clear() {
        m_elements.clear();
        m_vertices.clear();
    }

    const std::vector<MeshIO::IOElement> &elements() const { return m_elements; }
    const std::vector<MeshIO::IOVertex>  &vertices() const { return m_vertices; }

protected:
    std::vector<MeshIO::IOElement> m_elements;
    std::vector<MeshIO::IOVertex>  m_vertices;
};

template<size_t N>
class Inflator;

// 2D inflator is Luigi's WireInflator2D.
template<>
class Inflator<2> : public InflatorBase<Inflator<2>> {
public:
    typedef InflatorBase<Inflator<2>> Base;
    // TODO: make this piecewise linear instead of piecewise constant.
    typedef std::vector<Interpolant<Real, 1, 0>> NormalShapeVelocity;

    Inflator(const std::string &wireMeshPath)
        : m_inflator(wireMeshPath)
    {
        m_paramOp = m_inflator.patternGenerator().getParameterOperations();
    }

    void setMaxElementVolume(Real maxElementVol) { m_tparams.max_area = maxElementVol; }

    size_t numParameters() const { return m_inflator.patternGenerator().numberOfParameters(); }

    ParameterType parameterType(size_t p) const {
        ParameterType type;
        switch (m_paramOp.at(p).type) {
            case ParameterOperation::Radius:
                type = ParameterType::Thickness;
                break;
            case ParameterOperation::Translation:
                type = ParameterType::Offset;
                break;
            default: assert(false);
        }
        return type;
    }

    void inflate(const std::vector<Real> &params) {
        assert(params.size() == numParameters());
        CellParameters p_params = m_inflator.createParameters();
        for (size_t i = 0; i < params.size(); ++i)
            p_params.parameter(i) = params[i];

        if (!m_inflator.patternGenerator().parametersValid(p_params))
            throw runtime_error("Invalid parameters specified.");
        WireInflator2D::OutMeshType inflatedMesh;
        m_inflator.generatePattern(p_params, m_tparams, inflatedMesh);

        // Convert to MeshIO format
        clear();
        for (const auto &p : inflatedMesh.nodes)
            m_vertices.emplace_back(p[0], p[1], 0);
        for (const auto &e : inflatedMesh.elements)
            m_elements.emplace_back(e[0], e[1], e[2]);
        m_edge_fields = inflatedMesh.edge_fields;
    }

    // Needs access to mesh data structure to extract boundary fields
    // efficiently and dot velocities with normals
    template<class _FEMMesh>
    std::vector<NormalShapeVelocity>
    computeShapeNormalVelocities(const _FEMMesh &mesh) const {
        size_t numBE = mesh.numBoundaryElements();
        size_t nParams = numParameters();
        std::vector<NormalShapeVelocity> vn_p(nParams, NormalShapeVelocity(numBE));
        for (size_t bei = 0; bei < numBE; ++bei) {
            auto be = mesh.boundaryElement(bei);
            auto edge = make_pair(be.tip(). volumeVertex().index(),
                                  be.tail().volumeVertex().index());
            const auto &field = m_edge_fields.at(edge);
            assert(field.size() == nParams);
            // TODO: make a linear scalar field--velocity should be reported as
            // a per-vertex vector field that we dot with each boundary element
            // normal.
            for (size_t p = 0; p < nParams; ++p)
                vn_p[p][bei][0] = field[p];
        }

        return vn_p;
    }

private:
    // mutable because WireInflator2D doesn't have well-behaved constness
    mutable WireInflator2D m_inflator;
    TessellationParameters m_tparams;
    std::vector<ParameterOperation> m_paramOp;
    WireInflator2D::OutMeshType::EdgeFields m_edge_fields;
};

// 3D inflator is James' PeriodicExploration
template<>
class Inflator<3> : public InflatorBase<Inflator<3>> {
public:
    typedef InflatorBase<Inflator<3>> Base;
    // Piecewise linear normal shape velocity.
    typedef std::vector<Interpolant<Real, 2, 1>> NormalShapeVelocity;

    Inflator(const std::string &wireMeshPath)
        : m_inflator(wireMeshPath, 5.0, 0.5 * sqrt(2.0)) {
        m_inflator.with_all_parameters();
    }

    void setMaxElementVolume(Real maxElementVol) { m_maxElementVol = maxElementVol; }

    size_t numParameters() const { return m_inflator.get_num_dofs(); }

    ParameterType parameterType(size_t p) const {
        return m_inflator.is_thickness_dof(p) ? ParameterType::Thickness
                                              : ParameterType::Offset;
    }

    void inflate(const std::vector<Real> &params) {
        VectorF paramVector(params.size());
        for (size_t i = 0; i < params.size(); ++i) paramVector(i) = params[i];
        m_inflator.set_dofs(paramVector);
        try {
            m_inflator.periodic_inflate();
            m_inflator.run_tetgen(m_maxElementVol);
        }
        catch (...) {
            std::cerr << "Exception while inflating parameters:" << std::endl;
            for (size_t i = 0; i < params.size(); ++i) std::cerr << params[i] << "\t";
            std::cerr << std::endl;
            throw;
        }

        // Convert to MeshIO format.
        clear();
        MatrixIr els = m_inflator.get_voxels();
        for (size_t i = 0; i < size_t(els.rows()); ++i)
            m_elements.emplace_back(els(i, 0), els(i, 1), els(i, 2), els(i, 3));
        MatrixFr verts = m_inflator.get_vertices();
        for (size_t i = 0; i < size_t(verts.rows()); ++i)
            m_vertices.emplace_back(Point3D(verts.row(i)));
    }

    // Needs access to mesh data structure to dot velocities with normals
    template<class _FEMMesh>
    std::vector<NormalShapeVelocity>
    computeShapeNormalVelocities(const _FEMMesh &mesh) const {
        std::vector<MatrixFr> vertexVelocities = m_inflator.get_shape_velocities();
        size_t numBE = mesh.numBoundaryElements();
        size_t nParams = numParameters();

        std::vector<NormalShapeVelocity> vn_p(nParams, NormalShapeVelocity(numBE));
        for (size_t p = 0; p < nParams; ++p) {
            NormalShapeVelocity &vn = vn_p[p];
            const MatrixFr &vVel = vertexVelocities[p];
            for (size_t bei = 0; bei < numBE; ++bei) {
                auto be = mesh.boundaryElement(bei);
                assert(be.numVertices() == vn[bei].size());
                // Linearly interpolate normal velocity at each vertex.
                for (size_t n = 0; n < be.numVertices(); ++n)
                    vn[bei][n] = vVel.row(be.vertex(n).volumeVertex().index()).dot(be->normal());
            }
        }

        return vn_p;
    }

private:
    Real m_maxElementVol = 0.0;
    PeriodicExploration m_inflator;
};


#endif /* end of include guard: INFLATOR_HH */

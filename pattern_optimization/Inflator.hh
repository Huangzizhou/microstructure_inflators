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
#include <stdexcept>
#include <string>

namespace {
    VectorF ToVectorF(const std::vector<Real> &vec) {
        VectorF result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) result(i) = vec[i];
        return result;
    }
    std::vector<Real> FromVectorF(const VectorF &vec) {
        std::vector<Real> result(vec.rows());
        for (int i = 0; i < vec.rows(); ++i) result[i] = vec(i);
        return result;
    }
}

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
        setMaxElementVolume(0.0001);
    }

    void setMaxElementVolume(Real maxElementVol) { m_tparams.max_area = maxElementVol; }
    void configureSubdivision(const std::string &algorithm, size_t levels) {
        throw std::runtime_error("Subdivision not supported in 2D");
    }

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

    void setDoFOutputPrefix(const std::string &pathPrefix) {
        throw std::runtime_error("Writing pattern DoFs unsupported in 2D");
    }

    void writePatternDoFs(const std::string &path,
                          const std::vector<Real> &params) {
        throw std::runtime_error("Writing pattern DoFs unsupported in 2D");
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

    Inflator(const std::string &wireMeshPath,
             Real cell_size = 5.0, Real default_thickness = 0.5 * sqrt(2),
             bool isotropic_params = false)
        : m_inflator(wireMeshPath, cell_size, default_thickness) {
        if (isotropic_params) m_inflator.with_all_isotropic_parameters();
        else                  m_inflator.with_all_parameters();
        setMaxElementVolume(0.0);
    }

    void setMaxElementVolume(Real maxElementVol) { m_maxElementVol = maxElementVol; }
    void configureSubdivision(const std::string &algorithm, size_t levels) {
        m_inflator.with_refinement(algorithm, levels);
    }

    size_t numParameters() const { return m_inflator.get_num_dofs(); }

    ParameterType parameterType(size_t p) const {
        return m_inflator.is_thickness_dof(p) ? ParameterType::Thickness
                                              : ParameterType::Offset;
    }

    void inflate(const std::vector<Real> &params) {
        m_inflator.set_dofs(ToVectorF(params));
        try { m_inflate_dofs(); }
        catch (...) {
            std::cerr << "Exception while inflating parameters" << std::endl;
            for (size_t i = 0; i < params.size(); ++i) std::cerr << params[i] << "\t";
            std::cerr << std::endl;
            throw;
        }
    }

    void inflate(const std::string &dofFile) {
        m_inflator.load_dofs(dofFile);
        m_inflate_dofs();
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

    // Configure automatic logging of every set of inflated DoF parameters.
    // If pathPrefix is nonempty, a dof file will be written at
    // pathPrefix_$inflationNumber.dof
    void setDoFOutputPrefix(const std::string &pathPrefix) {
        m_dofOutPathPrefix = pathPrefix;
    }

    // Note, overwrites the dofs in m_inflator
    void loadPatternDoFs(const std::string &path, std::vector<Real> &params) {
        m_inflator.load_dofs(path);
        params = FromVectorF(m_inflator.get_dofs());
    }

    // Write parameters in James' DoF format.
    void writePatternDoFs(const std::string &path,
                          const std::vector<Real> &params) {
        m_inflator.set_dofs(ToVectorF(params));
        m_inflator.save_dofs(path);
    }

private:
    std::string m_dofOutPathPrefix;
    size_t m_inflationCount = 0;

    Real m_maxElementVol;
    PeriodicExploration m_inflator;

    // Inflate the DoFs already stored in the inflator.
    void m_inflate_dofs() {
        if (m_dofOutPathPrefix != "") {
            m_inflator.save_dofs(m_dofOutPathPrefix + "_" +
                    std::to_string(m_inflationCount) + ".dof");
        }

        m_inflator.periodic_inflate();
        m_inflator.run_tetgen(m_maxElementVol);

        ++m_inflationCount;

        // Convert to MeshIO format.
        clear();
        MatrixIr els = m_inflator.get_voxels();
        for (size_t i = 0; i < size_t(els.rows()); ++i)
            m_elements.emplace_back(els(i, 0), els(i, 1), els(i, 2), els(i, 3));
        MatrixFr verts = m_inflator.get_vertices();
        for (size_t i = 0; i < size_t(verts.rows()); ++i)
            m_vertices.emplace_back(Point3D(verts.row(i)));
    }
};


#endif /* end of include guard: INFLATOR_HH */

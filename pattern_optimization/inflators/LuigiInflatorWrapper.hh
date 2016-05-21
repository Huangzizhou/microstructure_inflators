#ifndef LUIGIINFLATORWRAPPER_HH
#define LUIGIINFLATORWRAPPER_HH

#include "../InflatorBase.hh"
// Must always include Functions.hh before WireInflator2D.h since vcg polutes global
// namespace
#include <Functions.hh> 
#include <memory>
#include "../../Luigi/wireinflator2D/src"

class LuigiInflatorWrapper : public InflatorBase<2> {
public:
    using NSV = NormalShapeVelocity<2>;

    Inflator(const std::string &wireMeshPath)
        : m_inflator(WireInflator2D::construct<WireMesh2D>(wireMeshPath))
    {
        m_paramOp = m_inflator->getParameterOperations();
        setMaxElementVolume(0.0001);
    }

    // MHS JUL14, 2015: 
    // A new constructor that takes in "const int symmetryMode" 
    // as its last parameter to pass to WireMesh2DMorteza 
    Inflator(const std::string &wireMeshPath, const int symmetryMode) {
        // symmetryMode < 0 reverts to Luigi's symmetry parameters (WireMesh2D)
        if (symmetryMode < 0) m_inflator = WireInflator2D::construct<WireMesh2D       >(wireMeshPath);
        else                  m_inflator = WireInflator2D::construct<WireMesh2DMorteza>(wireMeshPath, symmetryMode);
        m_paramOp = m_inflator->getParameterOperations();
        setMaxElementVolume(0.0001);
    }

    Inflator(const std::string &wireMeshPath,
             Real /* cell_size */, Real /* default_thickness = 0.5 * sqrt(2) */,
             bool /* isotropic_params = false */, bool /* vertex_thickness = false */)
        : m_inflator(WireInflator2D::construct<WireMesh2D>(wireMeshPath)) {
            throw std::runtime_error("2D inflator is not yet configurable");
    }

    void setMaxElementVolume(Real maxElementVol) { m_tparams.max_area = maxElementVol; }
    Real getMaxElementVolume() const { return m_tparams.max_area; }
    void configureSubdivision(const std::string &/* algorithm */, size_t /* levels */) {
        throw std::runtime_error("Subdivision not supported in 2D");
    }

    void setReflectiveInflator(bool use) { if (use) throw std::runtime_error("Luigi's inflator is not reflective"); }
    size_t numParameters() const { return m_inflator->numberOfParameters(); }

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
        CellParameters p_params = m_inflator->createParameters();
        for (size_t i = 0; i < params.size(); ++i)
            p_params.parameter(i) = params[i];

        if (!m_inflator->parametersValid(p_params))
            throw runtime_error("Invalid parameters specified.");
        WireInflator2D::OutMeshType inflatedMesh;
        m_inflator->generatePattern(p_params, m_tparams, inflatedMesh);

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
    std::vector<NSV>
    computeShapeNormalVelocities(const _FEMMesh &mesh) const {
        size_t numBE = mesh.numBoundaryElements();
        size_t nParams = numParameters();
        std::vector<NSV> vn_p(nParams, NSV(numBE));
        for (auto be : mesh.boundaryElements()) {
            auto edge = make_pair(be.tip(). volumeVertex().index(),
                                  be.tail().volumeVertex().index());
            const auto &field = m_edge_fields.at(edge);
            assert(field.size() == nParams);
            // TODO: make the 2D inflator output a true discontinuous piecwise
            // linear velocity field instead of its piecewise constant
            // approximation.
            for (size_t p = 0; p < nParams; ++p) {
                vn_p[p][be.index()][0] = field[p];
                vn_p[p][be.index()][1] = field[p];
            }
        }

        return vn_p;
    }

    // 2D is always printable.
    bool isPrintable(const std::vector<Real> &/* params */) const { return true; }

private:
    // mutable because WireInflator2D doesn't have well-behaved constness
    mutable std::shared_ptr<WireInflator2D> m_inflator;
    TessellationParameters m_tparams;
    std::vector<ParameterOperation> m_paramOp;
    WireInflator2D::OutMeshType::EdgeFields m_edge_fields;
};

#endif /* end of include guard: LUIGIINFLATORWRAPPER_HH */

#ifndef ISOINFLATORWRAPPER_HH
#define ISOINFLATORWRAPPER_HH

#include "../../isosurface_inflator/IsosurfaceInflator.hh"
#include "../../isosurface_inflator/IsosurfaceInflatorConfig.hh"
#include "../InflatorBase.hh"

template<size_t N>
class IsoinflatorWrapper : public InflatorBase<N> {
public:
    IsoinflatorWrapper(const std::string &wireMeshPath,
             Real /* cell_size */ = 5.0, Real /* default_thickness */ = 0.5 * sqrt(2),
             bool isotropic_params = false, bool vertex_thickness = false)
        : m_inflator(
                (N == 2) ? "2D_orthotropic" :
                (isotropic_params ? "cubic" : "orthotropic"),
                vertex_thickness, wireMeshPath)
    {
        // IsosurfaceInflatorConfig::get().inflationGraphPath = "inflgraph.wire";
    }

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
    virtual void inflate(const std::vector<Real> &params) {
        try {
            m_inflator.inflate(params);
            // Could avoid duplicate storage for this inflator...
            this->m_vertices = m_inflator.vertices();
            this->m_elements = m_inflator.elements();
        }
        catch (...) {
            std::cerr << setprecision(20);
            std::cerr << "Exception while inflating parameters" << std::endl;
            for (size_t i = 0; i < params.size(); ++i) std::cerr << params[i] << "\t";
            std::cerr << std::endl;
            throw;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    // Get a per-vertex perturbation vector field induced by changing
    // each parameter. Internal vertices get 0 velocity.
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const {
        std::vector<std::vector<Real>> nsv = m_inflator.normalShapeVelocities();
        std::vector<Point3D>             n = m_inflator.vertexNormals();

        const size_t np = numParameters();
        const size_t nv = n.size();
        assert(nv == nsv.at(0).size());

        std::vector<VectorField<Real, N>> result(np);
        for (size_t p = 0; p < np; ++p) {
            result[p].resizeDomain(nv);
            for (size_t vi = 0; vi < nv; ++vi) {
                result[p](vi)  = truncateFrom3D<VectorND<N>>(n[vi]);
                result[p](vi) *= nsv[p][vi];
                // assert(!std::isnan(result[p](vi)));
                // assert(!std::isnan(result[p](vi)));
            }
        }
        return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const    { return true; }
    virtual size_t numParameters() const { return m_inflator.numParams(); }
    virtual ParameterType parameterType(size_t p) const {
        if (m_inflator.isThicknessParam(p)) return ParameterType::Thickness;
        if (m_inflator. isPositionParam(p)) return ParameterType::Offset;
        if (m_inflator. isBlendingParam(p)) return ParameterType::Blending;
        throw std::runtime_error("Unknown parameter type for param " + std::to_string(p));
    }
    virtual bool isPrintable(const std::vector<Real> &params) const { return m_inflator.isPrintable(params); }

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) { m_inflator.meshingOptions().load(moptsPath); }
    virtual void setMaxElementVolume(Real maxElementVol) {
        if (N == 2) m_inflator.meshingOptions().maxArea = maxElementVol;
        else        m_inflator.meshingOptions().cellSize = maxElementVol;
    }
    virtual void setReflectiveInflator(bool use) { if (!use) throw std::runtime_error("IsosurfaceInflator is always reflective."); }

private:
    IsosurfaceInflator m_inflator;
};

#endif /* end of include guard: ISOINFLATORWRAPPER_HH */

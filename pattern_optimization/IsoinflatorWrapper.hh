#ifndef ISOINFLATORWRAPPER_HH
#define ISOINFLATORWRAPPER_HH

#include "../isosurface_inflator/IsosurfaceInflator.hh"

template<size_t N>
class IsoinflatorWrapper : public InflatorBase<IsoinflatorWrapper<N>> {
public:
    typedef InflatorBase<IsoinflatorWrapper<N>> Base;
    using NSV = NormalShapeVelocity<N>;

    IsoinflatorWrapper(const std::string &wireMeshPath,
             Real /* cell_size */ = 5.0, Real /* default_thickness */ = 0.5 * sqrt(2),
             bool isotropic_params = false, bool vertex_thickness = false)
        : m_inflator(
                (N == 2) ? "2D_orthotropic" :
                (isotropic_params ? "cubic" : "orthotropic"),
                vertex_thickness, wireMeshPath)
    { }

    void setMaxElementVolume(Real maxElementVol) {
        if (N == 2) m_inflator.meshingOptions().maxArea = maxElementVol;
        else        m_inflator.meshingOptions().cellSize = maxElementVol;
    }
    void configureSubdivision(const std::string &/* algorithm */, size_t levels) {
        if (levels != 0)
            throw std::runtime_error("IsosurfaceInflator doesn't support subdivision");
    }

    size_t numParameters() const { return m_inflator.numParams(); }

    ParameterType parameterType(size_t p) const {
        if (m_inflator.isThicknessParam(p)) return ParameterType::Thickness;
        if (m_inflator. isPositionParam(p)) return ParameterType::Offset;
        if (m_inflator. isBlendingParam(p)) return ParameterType::Blending;
        throw std::runtime_error("Unknown parameter type for param " + std::to_string(p));
    }

    void inflate(const std::vector<Real> &params) {
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

    void inflate(const std::string &/* dofFile */) {
        throw std::runtime_error("IsosurfaceInflator doesn't support DoF files.");
    }

    void setReflectiveInflator(bool use)      { if (!use) throw std::runtime_error("IsosurfaceInflator is always reflective."); }
    void setDumpSurfaceMesh(bool dump = true) { if (dump) throw std::runtime_error("IsosurfaceInflator does not do surface meshing."); }

    bool isPrintable(const std::vector<Real> &params) const {
        return m_inflator.isPrintable(params);
    }

    // NOTE: assumes periodic boundary conditions have already been applied
    // (for proper clearing of boundary velocity).
    template<class _FEMMesh>
    std::vector<NSV>
    computeShapeNormalVelocities(const _FEMMesh &mesh) const {
        // IsosurfaceInflator computes the velocity of each vertex as a scalar
        // normal velocity vnp, and a vertex normal.
        std::vector<std::vector<Real>> nsv_p = m_inflator.normalShapeVelocities();
        std::vector<Point3D>               n = m_inflator.vertexNormals();

        size_t numBE = mesh.numBoundaryElements();
        size_t nParams = numParameters();

        std::vector<NSV> vn_p(nParams, NSV(numBE));
        for (size_t p = 0; p < nParams; ++p) {
            NSV &vn = vn_p[p];
            const auto &nsv = nsv_p.at(p);
            for (auto be : mesh.boundaryElements()) {
                assert(be.numVertices() == vn[be.index()].size());
                if (be->isPeriodic) {
                    // Shape velocity should be exactly 0 on the periodic boundary.
                    // (If we don't enforce this, the nonzero shape velocity from
                    // (domega intersect dY) vertices will leak into the
                    // periodic boundary faces.
                    vn[be.index()] = 0;
                    continue;
                }
                // Interpolate the boundary element corner's normal velocities.
                for (size_t bvi = 0; bvi < be.numVertices(); ++bvi) {
                    size_t vi = be.vertex(bvi).volumeVertex().index();
                    // Note: the direct version empirically seems to be working better...
                    // It should be used by default.
                    if (PatternOptimization::Config::get().useSDNormalShapeVelocityDirectly)
                        vn[be.index()][bvi] = nsv.at(vi);
                    else {
                        // Get the normal velocity over each face by dotting with the face normal
                        // and linearly interpolating.
                        vn[be.index()][bvi] = nsv.at(vi) * truncateFrom3D<PointND<N>>(n.at(vi)).dot(be->normal());
                    }
                }
            }
        }

        return vn_p;
    }

    // Configure automatic logging of every set of inflated DoF parameters.
    // If pathPrefix is nonempty, a dof file will be written at
    // pathPrefix_$inflationNumber.dof
    void setDoFOutputPrefix(const std::string &/* pathPrefix */) {
        throw std::runtime_error("IsosurfaceInflator does not support DoF files.");
    }

    // Note, overwrites the dofs in m_inflator
    void loadPatternDoFs(const std::string &/* path */, std::vector<Real> &/* params */) {
        throw std::runtime_error("IsosurfaceInflator does not support DoF files.");
    }

    // Write parameters in James' DoF format.
    void writePatternDoFs(const std::string &/* path */,
                          const std::vector<Real> &/* params */) {
        throw std::runtime_error("IsosurfaceInflator does not support DoF files.");
    }

private:
    IsosurfaceInflator m_inflator;
};

#endif /* end of include guard: ISOINFLATORWRAPPER_HH */

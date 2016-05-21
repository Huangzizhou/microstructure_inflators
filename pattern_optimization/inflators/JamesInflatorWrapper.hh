#ifndef JAMESINFLATORWRAPPER_HH
#define JAMESINFLATORWRAPPER_HH

class JamesInaflatorWrapper : public InflatorBase<3> {
public:
    typedef InflatorBase<Inflator<3>> Base;
    using NSV = NormalShapeVelocity<3>;

    Inflator(const std::string &wireMeshPath,
             Real cell_size = 5.0, Real default_thickness = 0.5 * sqrt(2),
             bool isotropic_params = false, bool vertex_thickness = false)
        : m_inflator(wireMeshPath, cell_size, default_thickness) {
        ParameterCommon::TargetType thickness_type =
            vertex_thickness ? ParameterCommon::VERTEX : ParameterCommon::EDGE;
        if (isotropic_params) m_inflator.with_all_isotropic_parameters(thickness_type);
        else                  m_inflator.with_all_parameters(thickness_type);
        setMaxElementVolume(0.0);
    }

    void setMaxElementVolume(Real maxElementVol) { m_maxElementVol = maxElementVol; }
    Real getMaxElementVolume() const { return m_maxElementVol; }
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

    void setReflectiveInflator(bool use) { m_useReflectiveInflator = use; }
    void setDumpSurfaceMesh(bool dump = true) { m_dumpSurfaceMesh = dump; }

    // the printability check actually modifies the inflator's internal state,
    // so we can't mark this const.
    // Also, we must change the current dofs to run the printability test
    bool isPrintable(const std::vector<Real> &params) {
        m_inflator.set_dofs(ToVectorF(params));
        return m_inflator.is_printable();
    }

    // Needs access to mesh data structure to dot velocities with normals
    template<class _FEMMesh>
    std::vector<NSV>
    computeShapeNormalVelocities(const _FEMMesh &mesh) const {
        std::vector<MatrixFr> vertexVelocities = m_inflator.get_shape_velocities();
        size_t numBE = mesh.numBoundaryElements();
        size_t nParams = numParameters();

        std::vector<NSV> vn_p(nParams, NSV(numBE));
        for (size_t p = 0; p < nParams; ++p) {
            NSV &vn = vn_p[p];
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
    bool m_useReflectiveInflator = true;

    // Used for debugging when tetgen fails.
    bool m_dumpSurfaceMesh = false;

    // Inflate the DoFs already stored in the inflator.
    // (written to surface_debug.msh)
    void m_inflate_dofs() {
        if (m_dofOutPathPrefix != "") {
            m_inflator.save_dofs(m_dofOutPathPrefix + "_" +
                    std::to_string(m_inflationCount) + ".dof");
        }

        m_inflator.periodic_inflate(m_useReflectiveInflator);

        if (m_dumpSurfaceMesh) {
            // Debug surface mesh (for when tetgen fails)
            std::vector<MeshIO::IOElement> triangles;
            std::vector<MeshIO::IOVertex>  vertices;
            MatrixFr verts = m_inflator.get_vertices();
            MatrixIr facs = m_inflator.get_faces();
            for (size_t i = 0; i < size_t(facs.rows()); ++i)
                triangles.emplace_back(facs(i, 0), facs(i, 1), facs(i, 2));
            for (size_t i = 0; i < size_t(verts.rows()); ++i)
                vertices.emplace_back(Point3D(verts.row(i)));
            MeshIO::save("surface_debug.msh", vertices, triangles);
        }

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

#endif /* end of include guard: JAMESINFLATORWRAPPER_HH */

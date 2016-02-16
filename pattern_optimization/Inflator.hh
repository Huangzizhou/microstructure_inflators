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
// Must include Functions.hh before WireInflator2D.h since vcg polutes global
// namespace
#include <Functions.hh> 
#include <WireInflator2D.h>
#include <Wires/Interfaces/PeriodicExploration.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>
#include "rref.h"
#include "ParameterConstraint.hh"
#include "PatternOptimizationConfig.hh"
#include "InflatorBase.hh"

template<size_t N>
class Inflator;

#ifdef ISOSURFACE_INFLATOR
#include "IsoinflatorWrapper.hh"
// Use isosurface inflator as 2D inflator
template<>
class Inflator<2> : public IsoinflatorWrapper<2> {
public:
    using Base = IsoinflatorWrapper<2>;
    using Base::Base;
};
#else // !ISOSURFACE_INFLATOR
// Use Luigi's WireInflator2D as 2D inflator.
template<>
class Inflator<2> : public InflatorBase<Inflator<2>> {
public:
    typedef InflatorBase<Inflator<2>> Base;
    using NSV = NormalShapeVelocity<2>;

    Inflator(const std::string &wireMeshPath)
        : m_inflator(WireInflator2D::construct(wireMeshPath))
    {
        m_paramOp = m_inflator->getParameterOperations();
        setMaxElementVolume(0.0001);
    }

    // MHS JUL14, 2015: 
    // A new constructor that takes in "const int symmetryMode" 
    // as its last parameter to pass to WireMesh2DMorteza 
    Inflator(const std::string &wireMeshPath, const int symmetryMode)
        : m_inflator(WireInflator2D::construct(wireMeshPath, symmetryMode))
    {
        m_paramOp = m_inflator->getParameterOperations();
        setMaxElementVolume(0.0001);
    }


    Inflator(const std::string &wireMeshPath,
             Real /* cell_size */, Real /* default_thickness = 0.5 * sqrt(2) */,
             bool /* isotropic_params = false */, bool /* vertex_thickness = false */)
        : m_inflator(WireInflator2D::construct(wireMeshPath)) {
            throw std::runtime_error("2D inflator is not yet configurable");
    }

    void setMaxElementVolume(Real maxElementVol) { m_tparams.max_area = maxElementVol; }
    Real getMaxElementVolume() const { return m_tparams.max_area; }
    void configureSubdivision(const std::string &/* algorithm */, size_t /* levels */) {
        throw std::runtime_error("Subdivision not supported in 2D");
    }

    void setReflectiveInflator(bool /* use */) {
        throw std::runtime_error("Reflective inflator not supported in 2D");
    }

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

    void setDoFOutputPrefix(const std::string &/* pathPrefix */) {
        throw std::runtime_error("Writing pattern DoFs unsupported in 2D");
    }

    void writePatternDoFs(const std::string &/* path */,
                          const std::vector<Real> &/* params */) {
        throw std::runtime_error("Writing pattern DoFs unsupported in 2D");
    }

private:
    // mutable because WireInflator2D doesn't have well-behaved constness
    mutable std::shared_ptr<WireInflator2D> m_inflator;
    TessellationParameters m_tparams;
    std::vector<ParameterOperation> m_paramOp;
    WireInflator2D::OutMeshType::EdgeFields m_edge_fields;
};

#endif // !ISOSURFACE_INFLATOR

#ifdef ISOSURFACE_INFLATOR
#include "IsoinflatorWrapper.hh"
// Use isosurface inflator as 3D inflator
template<>
class Inflator<3> : public IsoinflatorWrapper<3> {
public:
    using Base = IsoinflatorWrapper<3>;
    using Base::Base;
};
#else // !ISOSURFACE_INFLATOR
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

// 3D inflator is James' PeriodicExploration
template<>
class Inflator<3> : public InflatorBase<Inflator<3>> {
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
#endif // !ISOSURFACE_INFLATOR

template <size_t N>
class ConstrainedInflator : public Inflator<N> {
public:
    typedef Inflator<N> Base;
    typedef typename Base::NSV NSV;

    template<typename... BaseArgs>
    ConstrainedInflator(const std::vector<string> &constraints,
                        BaseArgs... baseArgs)
        : Base(baseArgs...)
    {
        // Impose all equality constraints exactly
        // Solve for a set of independent variables, and express the
        // dependent variables in terms of these.
        size_t fullNumParams = Base::numParameters();
        size_t expectedCols = fullNumParams + 1;
        for (const auto &c : constraints) {
            ParameterConstraint pc(fullNumParams, c);
            if (!pc.isEqualityConstraint()) continue;
            m_augmentedConstraintSystem.emplace_back(pc.augmentedRow());
            if (m_augmentedConstraintSystem.back().size() != expectedCols)
                throw std::runtime_error("Invalid constraint row");
        }
        if (m_augmentedConstraintSystem.size() != constraints.size()) {
            std::cerr << "WARNING: Only equality constraints are currently implemented."
                      << std::endl;
        }

        // By construction, the augmented column will be made independent.
        // If this doesn't happen, the system is inconsistent.
        rref(m_augmentedConstraintSystem, m_depIdx, m_depRow);

        // Determine the independent variables (complement of dep)
        unsigned int currDepVar = 0;
        for (unsigned int i = 0; i < expectedCols; ++i) {
            if ((currDepVar < m_depIdx.size()) && (i == m_depIdx[currDepVar]))
                ++currDepVar;
            else m_indepIdx.push_back(i);
        }
        assert(currDepVar == m_depIdx.size());

        // The augmented column vector must have been made independent.
        if ((m_indepIdx.size() == 0) || (m_indepIdx.back() != expectedCols - 1))
            throw std::runtime_error("Couldn't make RHS independent");

        // The last column isn't actually a variable
        m_indepIdx.pop_back();
    }

    size_t numParameters() const { return m_indepIdx.size(); }

    void inflate(const std::vector<Real> &params) {
        Base::inflate(fullParametersForReduced(params));
    }

    void inflate(const std::string dofs) {
        if (m_depIdx.size() != 0) {
            std::cerr << "WARNING: inflating constrained pattern from DoFs that may not satisfy constraints!"
                << std::endl;
        }
        Base::inflate(dofs);
    }

    // Note: doesn't consider what dependent parameters p might control.
    ParameterType parameterType(size_t p) const {
        assert(p < m_indepIdx.size());
        return Base::parameterType(m_indepIdx[p]);
    }

    // Needs access to mesh data structure to dot velocities with normals
    template<class _FEMMesh>
    std::vector<NSV>
    computeShapeNormalVelocities(const _FEMMesh &mesh) const {
        std::vector<NSV> fullParamVelocity
            = Base::computeShapeNormalVelocities(mesh);

        size_t nParams = numParameters();
        std::vector<NSV> reducedParamVelocity;
        reducedParamVelocity.reserve(nParams);
        
        // Effectively apply transpose of change of variables matrix.
        for (size_t pIndep = 0; pIndep < nParams; ++pIndep) {
            // Get the independent variable's shape velocity.
            // (identity part of the change of variables matrix).
            size_t indepIdx = m_indepIdx[pIndep];
            reducedParamVelocity.push_back(fullParamVelocity.at(indepIdx));
            auto &vel = reducedParamVelocity.back();

            // Add in shape velocity for each variable dependent on pIndep.
            // (Chain rule)
            for (size_t pDep = 0; pDep < m_depIdx.size(); ++pDep) {
                const auto &row = m_augmentedConstraintSystem[m_depRow[pDep]];
                Real dPdep_dPindep = -row[indepIdx];
                if (std::abs(dPdep_dPindep) > 1e-10) {
                    const auto &depVel = fullParamVelocity.at(m_depIdx[pDep]);
                    for (size_t bei = 0; bei < depVel.size(); ++bei)
                        vel[bei] += dPdep_dPindep * depVel[bei];
                }
            }
        }

        return reducedParamVelocity;
    }

    // Forced non-const due to James' inflator interface.
    bool isPrintable(const std::vector<Real> &params) {
        return Base::isPrintable(fullParametersForReduced(params));
    }

    // Write parameters in James' DoF format.
    void writePatternDoFs(const std::string &path,
                          const std::vector<Real> &params) {
        Base::writePatternDoFs(path, fullParametersForReduced(params));
    }

    // Effectively apply the change of variables matrix.
    std::vector<Real> fullParametersForReduced(const std::vector<Real> &rparams) {
        std::vector<Real> fullParams(Base::numParameters());

        // Copy the independent params over
        // (identity part of change of variables matrix)
        assert(rparams.size() == m_indepIdx.size());
        for (size_t i = 0; i < m_indepIdx.size(); ++i) {
            assert(m_indepIdx[i] < fullParams.size());
            fullParams[m_indepIdx[i]] = rparams[i];
        }

        // Compute the dependent param values
        for (size_t i = 0; i < m_depIdx.size(); ++i) {
            Real param = 0;
            assert(m_depRow[i] < m_augmentedConstraintSystem.size());
            const auto &row = m_augmentedConstraintSystem[m_depRow[i]];

            for (size_t j = 0; j < m_indepIdx.size(); ++j)
                param -= rparams[j] * row[m_indepIdx[j]];

            // Include the constraint's RHS value.
            param += row.back();
            assert(m_depIdx[i] < fullParams.size());
            fullParams[m_depIdx[i]] = param;
        }

        return fullParams;
    }

private:
    // The full parameters corresponding to each reduced (independent) parameter
    std::vector<size_t> m_indepIdx;
    // The full parameters corresponding to the dependent parameters
    // (empty if no equality constraints are specified)
    std::vector<size_t> m_depIdx;
    // The row of m_augmentedConstraintSystem that holds each dependent
    // parameter's dependencies.
    std::vector<size_t> m_depRow;
    // Constraint system in reduced row echelon form (row major)
    // The last column holds the system RHS.
    std::vector<vector<Real>> m_augmentedConstraintSystem;
};

#endif /* end of include guard: INFLATOR_HH */

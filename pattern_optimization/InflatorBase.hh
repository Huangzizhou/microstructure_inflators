////////////////////////////////////////////////////////////////////////////////
// InflatorBase.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Abstract base class for inflators.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/10/2015 22:16:39
////////////////////////////////////////////////////////////////////////////////
#ifndef INFLATORBASE_HH
#define INFLATORBASE_HH

#include <vector>
#include <MeshIO.hh>
#include <Functions.hh> 
#include <stdexcept>

// Meta parameters are for EqualityConstrainedInflator
enum class ParameterType { Thickness, Offset, Blending, Meta };

template<size_t N>
class InflatorBase {
public:
    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
    virtual void inflate(const std::vector<Real> &params) = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Geometry access (dimension agnostic)
    ////////////////////////////////////////////////////////////////////////////
    const std::vector<MeshIO::IOElement> &elements() const { return m_elements; }
    const std::vector<MeshIO::IOVertex>  &vertices() const { return m_vertices; }
    void clear() { m_elements.clear(), m_vertices.clear(); }

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation/access
    // All inflators must implement per-volume-vertex parameter shape
    // velocities (volumeShapeVelocities). Then a per-boundary-vertex shape
    // velocity is extracted by ignoring the internal values.
    //
    // This is done because the inflators themselves might not know which are
    // boundary vertices (or might not have the same boundary vertex indexing)
    ////////////////////////////////////////////////////////////////////////////
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const = 0;

    // Extract boundary shape velocities from volume shape velocities.
    // Non-virtual! 
    template<class _FEMMesh>
    std::vector<VectorField<Real, N>> shapeVelocities(const _FEMMesh &mesh) const {
        auto vvels = volumeShapeVelocities();
        const size_t nvels = vvels.size();
        std::vector<VectorField<Real, N>> result(nvels, VectorField<Real, N>(mesh.numBoundaryVertices()));
        for (size_t i = 0; i < nvels; ++i) {
            result[i].clear();
            for (auto bv : mesh.boundaryVertices())
                result[i](bv.index()) = vvels[i](bv.volumeVertex().index());
        }
        return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const = 0;
    virtual size_t numParameters() const = 0;
    virtual ParameterType parameterType(size_t p) const = 0;
    virtual bool isPrintable(const std::vector<Real> &params) const = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) {
        throw std::runtime_error("This inflator doesn't support meshing options files");
    }

    virtual void setMaxElementVolume(Real maxElementVol) { throw std::runtime_error("Unimplemented"); }
    virtual Real getMaxElementVolume()                   { throw std::runtime_error("Unimplemented"); }
    virtual void setReflectiveInflator(bool use)         { throw std::runtime_error("Unimplemented"); }
    virtual void setDumpSurfaceMesh(bool dump = true)    { if (dump) throw std::runtime_error("This inflator does not support surface meshing."); }

    virtual void configureSubdivision(const std::string &/* algorithm */, size_t levels) {
        if (levels != 0)
            throw std::runtime_error("This inflator doesn't support subdivision");
    }

    virtual ~InflatorBase() { }

protected:
    std::vector<MeshIO::IOElement> m_elements;
    std::vector<MeshIO::IOVertex>  m_vertices;
};

#endif /* end of include guard: INFLATORBASE_HH */

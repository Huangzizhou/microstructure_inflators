////////////////////////////////////////////////////////////////////////////////
// Inflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Abstract base class for inflators.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/10/2015 22:16:39
////////////////////////////////////////////////////////////////////////////////
#ifndef INFLATOR_HH
#define INFLATOR_HH

#include <vector>
#include <MeshIO.hh>
#include <Functions.hh> 
#include <stdexcept>
#include <string>
#include <iostream> 
#include <iomanip>

#include <Fields.hh>

// Meta parameters are for EqualityConstrainedInflator
enum class ParameterType { Thickness, Offset, Blending, Meta };

inline std::string parameterTypeString(const ParameterType &type) {
    switch(type) {
        case ParameterType::Thickness:
            return "Thickness";
        case ParameterType::Offset:
            return "Offset";
        case ParameterType::Blending:
            return "Blending";
        case ParameterType::Meta:
            return "Meta";
        default:
            return "Invalid";
    }
}

// The dimension-independent part of the interface
class InflatorBase {
public:
    virtual size_t dimension() const = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
    virtual void inflate(const std::vector<Real> &params) = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Geometry access (dimension agnostic)
    ////////////////////////////////////////////////////////////////////////////
    virtual const std::vector<MeshIO::IOElement> &elements() const { return m_elements; }
    virtual const std::vector<MeshIO::IOVertex>  &vertices() const { return m_vertices; }
    virtual void clear() { m_elements.clear(), m_vertices.clear(); }

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const = 0;
    virtual size_t numParameters() const = 0;
    virtual ParameterType parameterType(size_t p) const = 0;
    // James' inflator's printability check modifies state :(
    virtual bool isPrintable(const std::vector<Real> &params) = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &/* moptsPath */) {
        throw std::runtime_error("This inflator doesn't support meshing options files");
    }

    virtual void setMaxElementVolume(Real /* maxElementVol */) { throw std::runtime_error("Unimplemented"); }
    virtual Real getMaxElementVolume() const                   { throw std::runtime_error("Unimplemented"); }
    virtual void setReflectiveInflator(bool /* use */)         { throw std::runtime_error("Unimplemented"); }
    virtual void setDumpSurfaceMesh(bool dump = true)          { if (dump) throw std::runtime_error("This inflator does not support surface meshing."); }

    virtual void configureSubdivision(const std::string &/* algorithm */, size_t levels) {
        if (levels != 0)
            throw std::runtime_error("This inflator doesn't support subdivision");
    }

    virtual ~InflatorBase() { }

protected:
    std::vector<MeshIO::IOElement> m_elements;
    std::vector<MeshIO::IOVertex>  m_vertices;

    void m_reportInflationException(const std::vector<Real> &params) const {
        std::cerr << std::setprecision(20);
        std::cerr << "Exception while inflating parameters" << std::endl;
        for (size_t i = 0; i < params.size(); ++i) std::cerr << params[i] << "\t";
        std::cerr << std::endl;
    }
};

// The dimension-dependent part of the interface
template<size_t N>
class Inflator : public InflatorBase {
public:
    virtual size_t dimension() const override { return N; }

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

    // Extract (true) boundary shape velocities from volume shape velocities.
    // Non-virtual! 
    template<class _FEMMesh>
    std::vector<VectorField<Real, N>> shapeVelocities(const _FEMMesh &mesh) const {
        auto vvels = volumeShapeVelocities();
        const size_t nvels = vvels.size();
        std::vector<VectorField<Real, N>> result(nvels, VectorField<Real, N>(mesh.numBoundaryVertices()));

        std::vector<bool> isTrueBoundaryVertex(mesh.numBoundaryVertices(), false);
        for (auto be : mesh.boundaryElements()) {
            if (be->isPeriodic) continue;
            for (auto bv : be.vertices())
                isTrueBoundaryVertex.at(bv.index()) = true;
        }

        for (size_t i = 0; i < nvels; ++i) {
            result[i].clear();
            for (auto bv : mesh.boundaryVertices()) {
                if (!isTrueBoundaryVertex.at(bv.index())) continue;
                result[i](bv.index()) = vvels[i](bv.volumeVertex().index());
            }
        }
        return result;
    }

    virtual ~Inflator() { }
};

#endif /* end of include guard: INFLATOR_HH */

////////////////////////////////////////////////////////////////////////////////
// VoxelsInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
// Inflator that, given a matrix with densities, produces a mesh
*/
////////////////////////////////////////////////////////////////////////////////

#ifndef RBFINFLATOR_HH
#define RBFINFLATOR_HH

#include "../Inflator.hh"
#include "IsoinflatorWrapper.hh"
#include <level_set/RBF.hh>
#include <isosurface_inflator/MeshingOptions.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <MeshFEM/FEMMesh.hh>
#include <MeshFEM/PeriodicBoundaryMatcher.hh>

#include <memory>
#include <utility>
#include <vector>
#include <string>
#include <limits>

class RBFInflator : public Inflator<2> {
public:
    using Mesh = FEMMesh<2, 1, VectorND<2>>;

    RBFInflator(Real epsilon, size_t dim);

    RBFInflator(std::string png_path, Real epsilon, size_t dim);

    ~RBFInflator() { }

    virtual MeshingOptions &meshingOptions() override { return m_meshingOptions; }

private:

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
    virtual void m_inflate(const std::vector<Real> &params) override;

public:

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    virtual std::vector<VectorField<Real, 2>> volumeShapeVelocities() const override;


    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override;
    virtual bool isPrintable(const std::vector<Real> &params) override { return false; }
    virtual ParameterType parameterType(size_t p) const override { return ParameterType::Custom1; }
    virtual std::vector<Real> defaultParameters() const override { return matToVec(m_coeffMatrix); }


    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override;
    virtual void setMaxElementVolume(Real maxElementVol) override;
    virtual Real getMaxElementVolume() const override;


    ////////////////////////////////////////////////////////////////////////////
    // RBFInflator-specific
    ////////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<Real>> vecToMat(const std::vector<Real> &vec, size_t nRows, size_t nCols) const;
    std::vector<Real> matToVec(const std::vector<std::vector<Real>> &mat) const;


    ////////////////////////////////////////////////////////////////////////////
    // Data members
    ////////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<Real>> m_coeffMatrix;
    MeshingOptions m_meshingOptions;
    size_t m_dim;
    Real m_epsilon;
    MidplaneMesher m_mesher;
    std::unique_ptr<Mesh> m_mesh;
    BBox<Point2D> m_bbox;
};
#endif /* end of include guard: RBFINFLATOR_HH */

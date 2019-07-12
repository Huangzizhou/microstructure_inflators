////////////////////////////////////////////////////////////////////////////////
// VoxelsInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
// Inflator that, given a matrix with densities, produces a mesh
*/
////////////////////////////////////////////////////////////////////////////////

#ifndef VOXELSINFLATOR_HH
#define VOXELSINFLATOR_HH

#include "../Inflator.hh"
#include "IsoinflatorWrapper.hh"
#include <voxel_inflator/VoxelSDF.hh>
#include <isosurface_inflator/MeshingOptions.hh>
#include <isosurface_inflator/MidplaneMesher.hh>

#include <memory>
#include <utility>
#include <vector>
#include <string>
#include <limits>

class VoxelsInflator : public Inflator<2> {
public:
    VoxelsInflator(const std::vector<std::vector<Real>> &densityMatrix);
    ~VoxelsInflator() { }

private:

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
    virtual void m_inflate(const std::vector<Real> &params) override;


    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    virtual std::vector<VectorField<Real, 2>> volumeShapeVelocities() const override;


    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override;

    virtual ParameterType parameterType(size_t p) const override { return ParameterType::Custom1; }


    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override;

    virtual void setMaxElementVolume(Real maxElementVol) override;

    virtual Real getMaxElementVolume() const override;


    ////////////////////////////////////////////////////////////////////////////
    // Data members
    ////////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<Real>> m_densityMatrix;
    MeshingOptions m_meshingOptions;
    size_t m_nCols, m_nRows;
    MidplaneMesher m_mesher;
};
#endif /* end of include guard: VOXELSINFLATOR_HH */

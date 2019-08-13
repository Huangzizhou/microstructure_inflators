////////////////////////////////////////////////////////////////////////////////
// RBFOrthoInflator.hh
////////////////////////////////////////////////////////////////////////////////

#ifndef RBFORTHOINFLATOR_HH
#define RBFORTHOINFLATOR_HH

#include "../Inflator.hh"
#include "RBFInflator.hh"
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

class RBFOrthoInflator : public Inflator<2> {
public:
    using Mesh = FEMMesh<2, 1, VectorND<2>>;

    RBFOrthoInflator(Real epsilon, size_t dim);

    RBFOrthoInflator(std::string png_path, Real epsilon, size_t dim);

    ~RBFOrthoInflator() { }

    virtual MeshingOptions &meshingOptions() override { return m_rbfInflator.meshingOptions(); }

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
    virtual std::vector<Real> defaultParameters() const override {
        std::vector<Real> result;
        std::vector<Real> originalParams = m_rbfInflator.defaultParameters();

        for (unsigned i=0; i < numParameters(); i++) {
            result.push_back(originalParams[reducedParamToAll(i)[0]]);
        }

        return result;
    }


    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override;
    virtual void setMaxElementVolume(Real maxElementVol) override;
    virtual Real getMaxElementVolume() const override;
    virtual void setReflectiveInflator(bool /*use*/)  override { }


    ////////////////////////////////////////////////////////////////////////////
    // RBFInflator-specific
    ////////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<Real>> vecToMat(const std::vector<Real> &vec, size_t nRows, size_t nCols) const;
    std::vector<Real> matToVec(const std::vector<std::vector<Real>> &mat) const;
    std::vector<size_t> reducedParamToAll(size_t param) const;
    size_t originalToReducedParam(size_t originalParam) const;
    void savePng(const std::vector<Real> &reducedParams, std::string png_path) const;


    ////////////////////////////////////////////////////////////////////////////
    // Data members
    ////////////////////////////////////////////////////////////////////////////
    size_t m_dim;
    RBFInflator m_rbfInflator;
};
#endif /* end of include guard: RBFORTHOINFLATOR_HH */

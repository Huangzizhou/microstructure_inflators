//
// Created by Davi Colli Tozoni on 3/26/18.
//

#ifndef MICROSTRUCTURES_CONSTRAINEDISOINFLATOR_H
#define MICROSTRUCTURES_CONSTRAINEDISOINFLATOR_H

#include "../Inflator.hh"
#include "IsoinflatorWrapper.hh"

#include <memory>
#include <utility>
#include <vector>
#include <string>
#include <limits>

template<size_t N>
class ConstrainedIsoinflator : public Inflator<N> {
public:
    ConstrainedIsoinflator(const std::string &wireMeshPath, const std::string &symmetryType, bool vertex_thickness, const std::vector<bool> &paramsMask, const std::vector<double> &originalParams, size_t inflationGraphRadius = 2);

    ~ConstrainedIsoinflator() { }

    ////////////////////////////////////////////////////////////////////////////
    // Geometry access (dimension agnostic)
    ////////////////////////////////////////////////////////////////////////////
    virtual const std::vector<MeshIO::IOElement> &elements() const override { return m_infl->elements(); }
    virtual const std::vector<MeshIO::IOVertex>  &vertices() const override { return m_infl->vertices(); }
    virtual void clear() override { m_infl->clear(); }


    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override;


    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    // Translate shape velocities from original to constrained parameters
    // (Basically, calls isoinflator function and filters out the fixed parameters.)
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const override;


    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override;
    virtual ParameterType parameterType(size_t p) const override;
    virtual bool isPrintable(const std::vector<Real> &params) override;


    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override { m_infl->loadMeshingOptions(moptsPath); }
    virtual void setMaxElementVolume(Real maxElementVol) override { m_infl->setMaxElementVolume(maxElementVol); }
    virtual Real getMaxElementVolume() const             override { return m_infl->getMaxElementVolume(); }
    virtual void setReflectiveInflator(bool use)         override { m_infl->setReflectiveInflator(use); }
    virtual void setDumpSurfaceMesh(bool dump = true)    override { m_infl->setDumpSurfaceMesh(dump); }
    virtual void configureSubdivision(const std::string &algorithm, size_t levels) override { m_infl->configureSubdivision(algorithm, levels); }
    virtual void setOrthoBaseCell(bool ortho)            override { m_infl->setOrthoBaseCell(ortho); }
    virtual Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> selfSupportingConstraints(const std::vector<double> &params) const override;


    ////////////////////////////////////////////////////////////////////////////
    // ConstrainedIsoinflator-specific
    ////////////////////////////////////////////////////////////////////////////
    // Effectively apply the change of variables
    template<typename T>
    std::vector<T> constrainedToFullParameters(const std::vector<T> &constrainedParameters) const;

    ////////////////////////////////////////////////////////////////////////////
    // Data members
    ////////////////////////////////////////////////////////////////////////////
    std::unique_ptr<IsoinflatorWrapper<N>> m_infl;
    std::vector<double> m_originalParams;
    std::vector<int> m_filteredParamsToOriginal;
    size_t m_numFilteredParams;
};

#endif //MICROSTRUCTURES_CONSTRAINEDISOINFLATOR_H

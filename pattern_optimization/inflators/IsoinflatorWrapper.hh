#ifndef ISOINFLATORWRAPPER_HH
#define ISOINFLATORWRAPPER_HH

#include "../Inflator.hh"

#include <memory>
#include <vector>
#include <string>
#include <stdexcept>

#include <Fields.hh>

class IsosurfaceInflator;

template<size_t N>
class IsoinflatorWrapper : public Inflator<N> {
public:
    IsoinflatorWrapper(const std::string &wireMeshPath,
                       bool isotropic_params, bool vertex_thickness);

    ////////////////////////////////////////////////////////////////////////////
    // Override geometry access
    ////////////////////////////////////////////////////////////////////////////
    virtual const std::vector<MeshIO::IOElement> &elements() const override;
    virtual const std::vector<MeshIO::IOVertex>  &vertices() const override;
    virtual void clear() override;

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
    virtual void inflate(const std::vector<Real> &params) override;

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    // Get a per-vertex perturbation vector field induced by changing
    // each parameter. Internal vertices get 0 velocity.
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const override;

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override    { return true; }
    virtual size_t numParameters() const override;
    virtual ParameterType parameterType(size_t p) const override;
    virtual bool isPrintable(const std::vector<Real> &params) override;

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override;
    virtual void setMaxElementVolume(Real maxElementVol) override;
    virtual void setReflectiveInflator(bool use) override { if (!use) throw std::runtime_error("IsosurfaceInflator is always reflective."); }

    // Out-of-line destructor needed due to incompletely-typed members.
    virtual ~IsoinflatorWrapper();
private:
    std::unique_ptr<IsosurfaceInflator> m_inflator;
};

#endif /* end of include guard: ISOINFLATORWRAPPER_HH */

////////////////////////////////////////////////////////////////////////////////
// ConstrainedInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Inflator that provides interface for fixing parameters to constants or
//    using an independent single value for multiple parameters
*/
////////////////////////////////////////////////////////////////////////////////

#ifndef MICROSTRUCTURES_CONSTRAINEDINFLATOR_H
#define MICROSTRUCTURES_CONSTRAINEDINFLATOR_H


#include "../Inflator.hh"
#include "BoundaryPerturbationInflator.hh"

#include <memory>
#include <utility>
#include <vector>
#include <string>
#include <limits>

template<size_t N>
class ConstrainedInflator : public Inflator<N> {
public:

    ConstrainedInflator(std::unique_ptr<Inflator<N>> inflator, const std::vector<bool> &paramsMask, const std::vector<int> &gluedParams = std::vector<int>()) : m_infl(std::move(inflator)) {
        m_originalParams = m_infl->defaultParameters();

        assert(m_infl->numParameters() == paramsMask.size());
        assert(m_infl->numParameters() == gluedParams.size());

        // Creates map from constrained parameters to original ones
        size_t originalIdx = 0;
        m_numFilteredParams = 0;
        for (; originalIdx < paramsMask.size(); originalIdx++) {
            if (!paramsMask[originalIdx]) {

                std::vector<size_t> correspondingOriginalVertices;
                if (gluedParams.size() > 0 && gluedParams[originalIdx] >= 0) {
                    // Check if glued parameters is already in filtered parameters list
                    if (gluedParams[originalIdx] < (int) originalIdx) {
                        continue;
                    }

                    correspondingOriginalVertices = {originalIdx, (size_t) gluedParams[originalIdx]};
                }
                else {
                    correspondingOriginalVertices = {originalIdx};
                }

                m_filteredParamsToOriginal.push_back(correspondingOriginalVertices);
                m_numFilteredParams++;
            }
        }
    }
    ~ConstrainedInflator() { }

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
    virtual void m_inflate(const std::vector<Real> &params) override {
        std::vector<Real> originalParams = constrainedToFullParameters(params);

        /*std::cout << "non-filtered p:";
        for (size_t i = 0; i < originalParams.size(); ++i) {
            if (i < 40)
                std::cout << "\t" << originalParams[i];
            else if (i == 40)
                std::cout << "..." << std::endl;
        }
        std::cout << std::endl;
        */

        m_infl->inflate(originalParams);

        //auto inflated_mesh = m_infl->mesh();
        //MSHFieldWriter perturbed_writer("inflated.msh", inflated_mesh, true);
    }
public:

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    // Translate shape velocities from original to constrained parameters
    // (Basically, calls isoinflator function and filters out the fixed parameters.)
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const override {
        std::vector<VectorField<Real, N>> result;
        std::vector<VectorField<Real, N>> originalVelocities = m_infl->volumeShapeVelocities();

        for (unsigned i=0; i < numParameters(); i++) {
            VectorField<Real, N> velocity = originalVelocities[m_filteredParamsToOriginal[i][0]];
            for (unsigned j=1; j < m_filteredParamsToOriginal[i].size(); j++) {
                velocity += originalVelocities[m_filteredParamsToOriginal[i][j]];
            }
            result.push_back(velocity);
        }

        return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override { return m_numFilteredParams; }
    virtual ParameterType parameterType(size_t p) const override { return m_infl->parameterType(m_filteredParamsToOriginal[p][0]); }
    virtual bool isPrintable(const std::vector<Real> &params) override { return m_infl->isPrintable(constrainedToFullParameters(params)); }
    virtual BBox<Vector3D> meshingCell() override { return m_infl->meshingCell(); }
    virtual std::vector<Real> defaultParameters() const override {
        std::vector<Real> originalParams = m_infl->defaultParameters();
        std::vector<Real> result;

        for (unsigned i=0; i < numParameters(); i++) {
            result.push_back(originalParams[m_filteredParamsToOriginal[i][0]]);
        }

        return result;
    }
    virtual ScalarField<Real> paramsFromBoundaryVField(const VectorField<Real, N> &bdryVField) const override {
        ScalarField<Real> originalParams = m_infl->paramsFromBoundaryVField(bdryVField);
        ScalarField<Real> result(numParameters());

        for (unsigned i=0; i < numParameters(); i++) {
            result[i] = originalParams[m_filteredParamsToOriginal[i][0]];
        }

        return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override { m_infl->loadMeshingOptions(moptsPath); }
    virtual MeshingOptions &meshingOptions() override { return m_infl->meshingOptions(); }
    virtual void setMaxElementVolume(Real maxElementVol) override { m_infl->setMaxElementVolume(maxElementVol); }
    virtual Real getMaxElementVolume() const             override { return m_infl->getMaxElementVolume(); }
    virtual void setReflectiveInflator(bool /*use*/)         override { } // does not seem to matter for this inflator
    virtual void setDumpSurfaceMesh(bool dump = true)    override { m_infl->setDumpSurfaceMesh(dump); }
    virtual void configureSubdivision(const std::string &algorithm, size_t levels) override { m_infl->configureSubdivision(algorithm, levels); }
    virtual void setOrthoBaseCell(bool ortho)            override { m_infl->setOrthoBaseCell(ortho); }
    virtual Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> selfSupportingConstraints(const std::vector<double> &/*params*/) const override {
        Eigen::MatrixXd C(0, numParameters() + 1);
        return C;
    };

    ////////////////////////////////////////////////////////////////////////////
    // ConstrainedIsoinflator-specific
    ////////////////////////////////////////////////////////////////////////////
    // Effectively apply the change of variables
    template<typename T>
    std::vector<T> constrainedToFullParameters(const std::vector<T> &constrainedParameters) const {
        std::vector<T> result(m_originalParams.size());

        // copy original vector
        for (unsigned i = 0; i < m_originalParams.size(); i++) {
            result[i] = m_originalParams[i];
        }

        // change only the constrained parameters using our map
        for (unsigned i = 0; i < numParameters(); i++) {
            for (unsigned j = 0; j < m_filteredParamsToOriginal[i].size(); j++) {
                result[m_filteredParamsToOriginal[i][j]] = constrainedParameters[i];
            }
        }

        return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Data members
    ////////////////////////////////////////////////////////////////////////////
    std::unique_ptr<Inflator<N>> m_infl;
    std::vector<double> m_originalParams;
    std::vector<std::vector<size_t>> m_filteredParamsToOriginal;
    size_t m_numFilteredParams;
};


#endif //MICROSTRUCTURES_CONSTRAINEDINFLATOR_H

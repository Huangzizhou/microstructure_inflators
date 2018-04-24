////////////////////////////////////////////////////////////////////////////////
// HexaPillarsInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Objective is to reduce the number of parameters to have a robust optimization.
//  In this inflator, four parameters shape the other original ones.
*/
//  Author:  Davi Colli Tozoni, dtozoni@nyu.edu
//  Company:  New York University
//  Created:  08/25/2016 13:50
////////////////////////////////////////////////////////////////////////////////

#ifndef HEXAPILLARSINFLATOR_HH
#define HEXAPILLARSINFLATOR_HH

#include "../Inflator.hh"
#include "IsoinflatorWrapper.hh"
#include "../../inflators/hex_inflator/hexlib.h"

#include <memory>
#include <utility>
#include <vector>
#include <string>
#include <limits>

class HexaPillarsInflator : public Inflator<2> {
public:
    HexaPillarsInflator(const std::vector<Real> &initial_params, double p2, char structure_type);
    ~HexaPillarsInflator() { }

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
    // Translate shape velocities from original to hexa pillars parameters (chain rule)
    // (Effectively apply the transpose of the change of variables matrix.)
    //TODO: discover how to implement it. Maybe automatic differentiation
    virtual std::vector<VectorField<Real, 2>> volumeShapeVelocities() const override;


    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override {
        unsigned num_parameters = m_structure_type == '+'? 2 : 7;
        cout << "Number of parameters: " << num_parameters << endl;
        return num_parameters;
    }

    // If only one full parameter is controlled, p's type is given by inflator.
    // Otherwise, p is a metaparameter
    virtual ParameterType parameterType(size_t p) const override {
        if (m_structure_type == '-')
            switch (p) {
                case 0: return ParameterType::Custom1;
                case 1: return ParameterType::Custom3;
                case 2: return ParameterType::Custom4;
                case 3: return ParameterType::Custom5;
                case 4: return ParameterType::Custom6;
                case 5: return ParameterType::Custom7;
                case 6: return ParameterType::Custom8;
                default: throw std::runtime_error("Custom parameter type in HexaPillarsInflator not recognizable.");
            }
        else {
            switch (p) {
                case 0: return ParameterType::Custom1;
                case 1: return ParameterType::Custom4;
                default: throw std::runtime_error("Custom parameter type in HexaPillarsInflator not recognizable.");
            }
        }
        return ParameterType::Meta;
    }

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

    virtual Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    selfSupportingConstraints(const std::vector<double> &params) const override;


    ////////////////////////////////////////////////////////////////////////////
    // HexaPillarsInflator-specific
    ////////////////////////////////////////////////////////////////////////////
    // Effectively apply the change of variables matrix.
    template<typename T>
    std::vector<T> hexaPillarsToFullParameters(const std::vector<T> &hexaPillarsParameters) const;
    void configureResolution(const std::vector<Real> &params);


    ////////////////////////////////////////////////////////////////////////////
    // Data members
    ////////////////////////////////////////////////////////////////////////////
    std::unique_ptr<IsoinflatorWrapper<2>> m_infl;
    Real m_p1, m_p2, m_p3, m_p4, m_p5, m_p6, m_p7, m_p8;
    char m_structure_type;
};
#endif /* end of include guard: HEXAPILLARSINFLATOR_HH */

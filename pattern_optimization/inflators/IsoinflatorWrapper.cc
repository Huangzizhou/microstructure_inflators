#include "IsoinflatorWrapper.hh"

#include "../../isosurface_inflator/IsosurfaceInflator.hh"
#include "../../isosurface_inflator/IsosurfaceInflatorConfig.hh"
#include <Future.hh>

#include <iostream>
#include <iomanip>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
// Override geometry access
////////////////////////////////////////////////////////////////////////////////
template<size_t N> const std::vector<MeshIO::IOElement> &IsoinflatorWrapper<N>::elements() const { return m_inflator->elements(); }
template<size_t N> const std::vector<MeshIO::IOVertex>  &IsoinflatorWrapper<N>::vertices() const { return m_inflator->vertices(); }
template<size_t N> void IsoinflatorWrapper<N>::clear() { m_inflator->clear(); }

////////////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
IsoinflatorWrapper<N>::IsoinflatorWrapper(const std::string &wireMeshPath,
                                          bool isotropic_params, bool vertex_thickness) {
    m_inflator = Future::make_unique<IsosurfaceInflator>(
                (N == 2) ? "2D_orthotropic" :
                (isotropic_params ? "cubic" : "orthotropic"),
                vertex_thickness, wireMeshPath);
    if ((N == 2) && isotropic_params)
        throw std::runtime_error("2D isosurface inflator currently only supports orthotropic parameters.");
    // IsosurfaceInflatorConfig::get().inflationGraphPath = "inflgraph.wire";
}

////////////////////////////////////////////////////////////////////////////////
// Inflation
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
void IsoinflatorWrapper<N>::m_inflate(const std::vector<Real> &params)
{
    m_inflator->inflate(params);
}

////////////////////////////////////////////////////////////////////////////////
// Shape velocity computation
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
std::vector<VectorField<Real, N>>
IsoinflatorWrapper<N>::volumeShapeVelocities() const
{
    std::vector<std::vector<Real>> nsv = m_inflator->normalShapeVelocities();
    std::vector<Point3D>             n = m_inflator->vertexNormals();

    const size_t np = numParameters();
    const size_t nv = this->vertices().size();
    assert(nv == n.size());
    assert(nv == nsv.at(0).size());

    std::vector<VectorField<Real, N>> result(np);
    for (size_t p = 0; p < np; ++p) {
        result[p].resizeDomain(nv);
        for (size_t vi = 0; vi < nv; ++vi) {
            result[p](vi)  = truncateFrom3D<VectorND<N>>(n[vi]);
            result[p](vi) *= nsv[p][vi];
            // assert(!std::isnan(result[p](vi)));
            // assert(!std::isnan(result[p](vi)));
        }
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Queries
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
size_t
IsoinflatorWrapper<N>::numParameters() const { return m_inflator->numParams(); }

template<size_t N>
ParameterType
IsoinflatorWrapper<N>::parameterType(size_t p) const
{
    if (m_inflator->isThicknessParam(p)) return ParameterType::Thickness;
    if (m_inflator-> isPositionParam(p)) return ParameterType::Offset;
    if (m_inflator-> isBlendingParam(p)) return ParameterType::Blending;
    throw std::runtime_error("Unknown parameter type for param " + std::to_string(p));
}

template<size_t N>
bool IsoinflatorWrapper<N>::isPrintable(const std::vector<Real> &params) { return m_inflator->isPrintable(params); }

////////////////////////////////////////////////////////////////////////////////
// Configuration
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
void IsoinflatorWrapper<N>::loadMeshingOptions(const std::string &moptsPath) { m_inflator->meshingOptions().load(moptsPath); }

template<size_t N>
void IsoinflatorWrapper<N>::setMaxElementVolume(Real maxElementVol)
{
    if (N == 2) m_inflator->meshingOptions().maxArea = maxElementVol;
    else        m_inflator->meshingOptions().cellSize = maxElementVol;
}

template<size_t N> void IsoinflatorWrapper<N>::setOrthoBaseCell(bool ortho) { m_inflator->setGenerateFullPeriodCell(!ortho); }
template<size_t N> BaseCellType IsoinflatorWrapper<N>::baseCellType() const { return m_inflator->baseCellType(); } 

////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
IsoinflatorWrapper<N>::~IsoinflatorWrapper() { }

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: 2D and 3D inflators.
////////////////////////////////////////////////////////////////////////////////
template class IsoinflatorWrapper<2>;
template class IsoinflatorWrapper<3>;

#include "BoxIntersectionMesher.hh"

#include <vector>
#include <utility>

#include "WireMesh.hh"
#include "PatternSignedDistance.hh"
#include "BoxIntersection1DFeatures.hh"

using namespace std;

template<class SignedDistanceFunction>
void BoxIntersectionMesher<SignedDistanceFunction>::
mesh(const SignedDistanceFunction &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements)
{
    vector<vector<pair<size_t, size_t>>> polygons;
    vertices.clear();
    boxIntersection1DFeatures(sdf, meshingOptions.marchingSquaresGridSize,
                              meshingOptions.marchingSquaresCoarsening,
                              vertices, polygons);
    elements.clear(), elements.reserve(polygons.size());
    for (const auto &p : polygons)
        for (const auto &e : p)
            elements.emplace_back(e.first, e.second);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template class BoxIntersectionMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>>>;
// Enable for slower builds...
// template class BoxIntersectionMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>>>;
// template class BoxIntersectionMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>>>;

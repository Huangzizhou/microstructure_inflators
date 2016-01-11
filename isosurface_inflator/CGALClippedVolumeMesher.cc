#include "CGALClippedVolumeMesher.hh"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <vector>

#include "BoxIntersection1DFeatures.hh"
#include "WireMesh.hh"
#include "PatternSignedDistance.hh"

// avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Domain
typedef K::FT FT;
typedef K::Point_3 Point;

template<class SignedDistanceFunction>
struct CGALClippedVolumeMesher<SignedDistanceFunction>::
ClippedSignedDistanceFunction {
    ClippedSignedDistanceFunction(const SignedDistanceFunction &sdf)
        : m_sdf(sdf), m_meshingBox(SignedDistanceFunction::boundingBox()) { }

    FT operator()(const Point &p) const {
        Point3<Real> pp(p[0], p[1], p[2]);
        return std::max(m_sdf.signedDistance(pp),
                        m_meshingBox.signedDistance(pp));
    }
private:
    const SignedDistanceFunction &m_sdf;
    SD::Primitives::Box<Real> m_meshingBox;
};

template<class SignedDistanceFunction>
void CGALClippedVolumeMesher<SignedDistanceFunction>::
mesh(const SignedDistanceFunction &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements)
{
    typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Implicit_mesh_domain_3<ClippedSignedDistanceFunction,K> > Mesh_domain;

    // Polyline
    typedef std::vector<Point>    Polyline_3;
    typedef std::list<Polyline_3> Polylines;

    // Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
    typedef typename CGAL::Mesh_triangulation_3<
            Mesh_domain,
            CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
            CGAL::Parallel_tag                        // Tag to activate parallelism
        >::type Tr;
#else
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, typename Mesh_domain::Corner_index,
                                                    typename Mesh_domain::Curve_segment_index> C3t3;
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

    ClippedSignedDistanceFunction cgal_sdfunc(sdf);
    std::list<std::vector<MeshIO::IOVertex>> polylinesMeshIO;
    boxIntersection1DFeatures(sdf, meshingOptions.marchingSquaresGridSize, polylinesMeshIO);
    // Check for near-intersection of polylines--this should only happen if
    // we failed to stitch up the boundary curves correctly.
    // This brute-force O(n^2) could easily be sped up...
    for (auto it1 = polylinesMeshIO.begin(); it1 != polylinesMeshIO.end(); ++it1) {
        for (auto it2 = polylinesMeshIO.begin(); it2 != polylinesMeshIO.end(); ++it2) {
            // Close points on the same curve are fine.
            if (it1 == it2) continue;
            for (const auto &v1 : *it1) {
                for (const auto &v2 : *it2) {
                    if (((v1.point - v2.point).norm() < 1e-8) &&
                            ((v1.point[0] != v2.point[0]) ||
                             (v1.point[1] != v2.point[1]) ||
                             (v1.point[2] != v2.point[2]))) {
                        throw std::runtime_error("Inexact intersection of polylines"); 
                    }
                }
            }
        }
    }

    // Convert to CGAL's polyline format
    Polylines polylines;
    for (const auto &l : polylinesMeshIO) {
        polylines.push_back(Polyline_3());
        Polyline_3 &line = polylines.back();
        line.reserve(l.size());
        for (size_t i = 0; i < l.size(); ++i)
            line.emplace_back(l[i][0], l[i][1], l[i][2]);
    }

    Point3d c;
    double r;
    sdf.boundingSphere(c, r);

    Mesh_domain domain(cgal_sdfunc,
            K::Sphere_3(Point(c[0], c[1], c[2]), r * r), meshingOptions.domainErrorBound);
    // std::cout << "Adding features..." << std::endl;
    domain.add_features(polylines.begin(), polylines.end());

    // Mesh criteria
    Mesh_criteria criteria(facet_angle            = meshingOptions.facetAngle,
                           facet_size             = meshingOptions.facetSize,
                           facet_distance         = meshingOptions.facetDistance,
                           cell_radius_edge_ratio = meshingOptions.cellRadiusEdgeRatio,
                           cell_size              = meshingOptions.cellSize,
                           edge_size              = meshingOptions.edgeSize);

    // Mesh generation
    // std::cout << "Making mesh..." << std::endl;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    // Access triangulation directly
    const Tr &tr = c3t3.triangulation();
    vertices.clear(), elements.clear();
    vertices.reserve(tr.number_of_vertices());
    std::map<typename C3t3::Vertex_handle, size_t> V;
    size_t i = 0;
    for (auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it) {
        V[it] = i++;
        Point p = it->point();
        vertices.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()),
                CGAL::to_double(p.z()));
    }

    elements.reserve(c3t3.number_of_cells_in_complex());
    for (auto it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it) {
        elements.emplace_back(V[it->vertex(0)],
                              V[it->vertex(1)],
                              V[it->vertex(2)],
                              V[it->vertex(3)]);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template class CGALClippedVolumeMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>>>;
// Enable for slower builds...
// template class CGALClippedVolumeMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>>>;
// template class CGALClippedVolumeMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>>>;
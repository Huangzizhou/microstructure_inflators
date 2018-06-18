#ifndef CONVEXHULL_HH
#define CONVEXHULL_HH

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

#include <cassert>
#include <sstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_items_with_id_3          Items;
typedef CGAL::Polyhedron_3<K, Items>              Polyhedron_3;
typedef K::Point_3                                Point_3;

template<class PointCollection>
void convexHull(const PointCollection &points,
                std::vector<MeshIO::IOVertex > &hullVertices,
                std::vector<MeshIO::IOElement> &hullElements) {
    std::vector<Point_3> cgal_pts;
    const size_t npts = points.size();
    cgal_pts.reserve(npts);
    for (const auto &pt : points)
        cgal_pts.emplace_back(pt[0], pt[1], pt[2]);

    Polyhedron_3 hull;
    CGAL::convex_hull_3(cgal_pts.begin(), cgal_pts.end(), hull);
    assert(hull.is_closed());
    assert(hull.is_pure_triangle());

    hullVertices.clear();
    hullVertices.reserve(hull.size_of_vertices());
    for (auto it = hull.vertices_begin(); it != hull.vertices_end(); ++it) {
        it->id() = hullVertices.size();
        const Point_3& p = it->point();
        hullVertices.emplace_back(p.x(), p.y(), p.z());
    }

    hullElements.clear();
    hullElements.reserve(hull.size_of_facets());
    for (auto it = hull.facets_begin(); it != hull.facets_end(); it++) {
        auto h_it = it->facet_begin();
        MeshIO::IOElement e;
        do {
            e.push_back(h_it->vertex()->id());
        } while (++h_it != it->facet_begin());
        assert(e.size() == 3);
        hullElements.push_back(e);
    }
}


#endif /* end of include guard: CONVEXHULL_HH */

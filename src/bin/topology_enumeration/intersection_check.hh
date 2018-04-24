#ifndef INTERSECTION_CHECK_HH
#define INTERSECTION_CHECK_HH

#include "../isosurface_inflator/Symmetry.hh"
#include <vector>
#include <utility>
#include <MeshIO.hh>

bool hasSelfIntersection(const std::vector<Point3d> &nodes,
                         const std::vector<std::pair<size_t, size_t>> &edges);

#endif /* end of include guard: INTERSECTION_CHECK_HH */

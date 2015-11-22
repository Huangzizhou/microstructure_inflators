////////////////////////////////////////////////////////////////////////////////
// CGALClippedVolumeMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      CGAL-based volume mesher for meshing a smooth signed distance function
//      intersected with an axis-aligned box; assumes that the only sharp
//      features come from the boolean intersection, and those sharp feature
//      curves are resolved using marching squares on the box faces.
//
//      The SignedDistanceFunction template parameter should be a class
//      provding the following functions:
//          static BBox<Point3<Real>> boundingBox()
//          template<R> signedDistance(Point3<R> p) const
//          template<R> boundingSphere(Point3<R> &c, R &r) const
//      boundingBox() determines the clipping region.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/29/2015 15:42:57
////////////////////////////////////////////////////////////////////////////////
#ifndef CGALCLIPPEDVOLUMEMESHER_HH
#define CGALCLIPPEDVOLUMEMESHER_HH

#include <vector>
#include <MeshIO.hh>
#include "MeshingOptions.hh"

template<class SignedDistanceFunction>
class CGALClippedVolumeMesher {
public:
    typedef typename SignedDistanceFunction::Real Real;
    CGALClippedVolumeMesher(const MeshingOptions &opts = MeshingOptions())
        : meshingOptions(opts) { }

    void mesh(const SignedDistanceFunction &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements);

    MeshingOptions meshingOptions;
private:
    struct ClippedSignedDistanceFunction;
};

#endif /* end of include guard: CGALCLIPPEDVOLUMEMESHER_HH */

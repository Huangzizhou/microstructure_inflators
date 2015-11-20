////////////////////////////////////////////////////////////////////////////////
// BoxIntersectionMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Mesher of the 1D intersection of the bounding box with the domain's
//      surface (intended for debugging).
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  11/12/2015 13:06:17
////////////////////////////////////////////////////////////////////////////////
#ifndef BOXINTERSECTIONMESHER_HH
#define BOXINTERSECTIONMESHER_HH

#include "MeshingOptions.hh"
#include <MeshIO.hh>

template<class SignedDistanceFunction>
class BoxIntersectionMesher {
public:
    typedef typename SignedDistanceFunction::Real Real;
    BoxIntersectionMesher(const MeshingOptions &opts = MeshingOptions())
        : meshingOptions(opts) { }

    void mesh(const SignedDistanceFunction &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements);

    MeshingOptions meshingOptions;
};

#endif /* end of include guard: BOXINTERSECTIONMESHER_HH */

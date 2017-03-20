////////////////////////////////////////////////////////////////////////////////
// MidplaneMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Mesh midplane slice through the 3D signed distance function. I.e. mesh
//      the 1D intersection of the 3D pattern surface with the z = 0 midplane.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  02/03/2016 11:41:48
////////////////////////////////////////////////////////////////////////////////
#ifndef MIDPLANEMESHER_HH
#define MIDPLANEMESHER_HH

#include "MeshingOptions.hh"
#include <MeshIO.hh>
#include "SignedDistanceRegion.hh"

class MidplaneMesher {
public:
    MidplaneMesher(const MeshingOptions &opts = MeshingOptions())
        : meshingOptions(opts) { }

    void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements);

    MeshingOptions meshingOptions;
};

#endif /* end of include guard: MIDPLANEMESHER_HH */

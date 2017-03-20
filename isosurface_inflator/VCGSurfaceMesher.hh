////////////////////////////////////////////////////////////////////////////////
// VCGSurfaceMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      VCG-based surface mesher (used for quick previews of the object).
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  11/11/2015 20:29:11
////////////////////////////////////////////////////////////////////////////////
#ifndef VCGSURFACEMESHER_HH
#define VCGSURFACEMESHER_HH

#include "MeshingOptions.hh"
#include <MeshIO.hh>

#include "SignedDistanceRegion.hh"

class VCGSurfaceMesher {
public:
    using Real = SignedDistanceRegion<3>::Real;
    VCGSurfaceMesher(const MeshingOptions &opts = MeshingOptions())
        : meshingOptions(opts) { }

    void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements);

    MeshingOptions meshingOptions;
};

#endif /* end of include guard: VCGSURFACEMESHER_HH */

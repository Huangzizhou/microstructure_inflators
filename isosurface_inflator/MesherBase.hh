////////////////////////////////////////////////////////////////////////////////
// MesherBase.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Abstract base class for all isosurface inflator meshers.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/06/2017 02:27:47
////////////////////////////////////////////////////////////////////////////////
#ifndef MESHERBASE_HH
#define MESHERBASE_HH

#include "MeshingOptions.hh"
#include <MeshIO.hh>
#include "SignedDistanceRegion.hh"

class MesherBase {
public:
    virtual void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements) = 0;

    MeshingOptions meshingOptions;
    bool periodic = false;
};

#endif /* end of include guard: MESHERBASE_HH */

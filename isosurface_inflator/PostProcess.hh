////////////////////////////////////////////////////////////////////////////////
// PostProcess.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Postprocess/clean up the geometry from the mesher and then compute
//  normals and shape velocities.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  08/23/2017 16:49:49
////////////////////////////////////////////////////////////////////////////////
#ifndef POSTPROCESS_HH
#define POSTPROCESS_HH

#include <vector>
#include <MeshIO.hh>
#include <Geometry.hh>
#include "IsosurfaceInflator.hh"

// Postprocess:
//    Snap to base cell and then reflect if necessary
//    Compute vertex normals and normal shape velocities
template<size_t N, class Point>
void postProcess(std::vector<MeshIO::IOVertex>  &vertices,
                 std::vector<MeshIO::IOElement> &elements,
                 std::vector<std::vector<Real>> &normalShapeVelocities,
                 std::vector<Point>             &vertexNormals,
                 const IsosurfaceInflator::Impl &inflator,
                 bool                      meshedFullPeriodCell,
                 bool                      requestFullPeriodCell,
                 const BBox<Point>         &meshCell,
                 const MeshingOptions      &opts);

#endif /* end of include guard: POSTPROCESS_HH */

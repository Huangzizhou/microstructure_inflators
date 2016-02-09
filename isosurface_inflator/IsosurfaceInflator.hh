////////////////////////////////////////////////////////////////////////////////
// IsosurfaceInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Interface to the isosurface inflator.
//      Everything is implemented in terms of opaque class
//      IsosurfaceInflatorImpl to speed up compilation of the code including us
//      (compiling anything depending on CGAL is slow!!).
//
//      Both 2D and 3D isosurface inflators are implemented (and can be selected
//      using the constructor's "type" argument, but the result is always
//      embedded in 3D.
//      
//      This inflator remembers the inflated geometry 
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/29/2015 12:17:26
////////////////////////////////////////////////////////////////////////////////
#ifndef ISOSURFACEINFLATOR_HH
#define ISOSURFACEINFLATOR_HH
#include <string>
#include <vector>
#include <MeshIO.hh>
#include "WireMesh.hh"
#include "MeshingOptions.hh"

class IsosurfaceInflator {
public:
    typedef Eigen::Matrix<Real, 3, 1> Point;
    enum class ParameterType { Thickness, Position, Blending };

    IsosurfaceInflator(const std::string &type, bool vertexThickenss,
                       const std::string &wireMeshPath);

    void inflate(const std::vector<Real> &params);
    std::vector<Real> defaultParameters() const;

    // Configure whether the full period cell is generated. E.g. in the case of
    // orthotropic symmetry, controls whether the representative eighth-cell is
    // reflected into the full cell.
    void setGenerateFullPeriodCell(bool onoff);
    bool shouldGenerateFullPeriodCell() const;

    // Access the inflation result
    // Vertex normals and shape velocity are per-vertex, taking value 0 on
    // non-boundary vertices.
    const std::vector<MeshIO::IOVertex>  &vertices() const;
    const std::vector<MeshIO::IOElement> &elements() const;
    const std::vector<std::vector<Real>> &normalShapeVelocities() const;
    const std::vector<Point>             &vertexNormals() const;

    size_t numParams() const;

    bool isThicknessParam(size_t p) const;

    // For debugging
    void disablePostprocess();
    void  enablePostprocess();

    MeshingOptions &meshingOptions();

    bool isPrintable(const std::vector<Real> &/* params */) const {
        // TODO: implement!
        return true;
    }

    ~IsosurfaceInflator();

public:
    class Impl;
    Impl *m_imp;
};

#endif /* end of include guard: ISOSURFACEINFLATOR_HH */

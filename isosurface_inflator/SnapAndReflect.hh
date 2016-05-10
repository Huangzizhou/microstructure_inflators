////////////////////////////////////////////////////////////////////////////////
// SnapAndReflect.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Snap representative base cell geometry to exactly fill [0, 1]^3, then
//      reflect into [-1, 1]^3 via the orthotropic symmetry planes.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/30/2015 14:37:00
////////////////////////////////////////////////////////////////////////////////
#ifndef SNAPANDREFLECT_HH
#define SNAPANDREFLECT_HH

#include <vector>
#include <ratio>
#include <MeshIO.hh>
#include "Isometries.hh"

inline bool isEq(Real a, Real b, Real tol = 0) {
    return std::abs(a - b) < tol;
}

// TODO: choose tolerance that works with both 2D and 3D inflator?
template<typename Vertex, typename TOL = std::ratio<2, long(1e4)>>
void snapVerticesToUnitCell(std::vector<Vertex> &vertices,
                            std::vector<std::vector<bool>> &onMinFace,
                            std::vector<std::vector<bool>> &onMaxFace) {
    static constexpr double tolerance = double(TOL::num) / double(TOL::den);
    // Snap vertices to [0, 1], determining min/max face membership
    onMinFace.assign(3, std::vector<bool>(vertices.size(), false));
    onMaxFace.assign(3, std::vector<bool>(vertices.size(), false));

    for (size_t vi = 0; vi < vertices.size(); ++vi) {
        auto &v = vertices[vi];
        for (size_t c = 0; c < 3; ++c) {
            if (isEq(v[c], 0, tolerance)) { v[c] = 0; onMinFace[c][vi] = true; }
            if (isEq(v[c], 1, tolerance)) { v[c] = 1; onMaxFace[c][vi] = true; }
        }
    }
}

// Generate full reflected cell in three steps: reflect along x axis,
// reflect along y axis, then reflect along z
//     +------+------+
//    /  reflect z  /|
//   +------+------+ |
//  /      /      /| |
// +------+------+ | +
// |      |      | |/|
// | refx | base | + |
// |      | (0)  |/| |
// +------+------+ | |
// |      |      | |/ 
// |  reflect y  | |
// |      |      |/
// +------+------+
// Vertices on the reflection planes must not be duplicated, and the
// onReflectionPlane arrays allow us to enforce this.
// onReflectionPlane[d] stores whether vertices lie on reflection plane d.
// We update these arrays with each reflection to keep them valid.
// If Dim = 2, only X and Y reflection are performed.
template<typename Vertex, typename Element>
void reflectXYZ(size_t Dim, // Dimensions to reflect in (length of [x, y, z] prefix)
                const std::vector<Vertex> &vertices,
                const std::vector<Element> &elements,
                std::vector<std::vector<bool>> onReflectionPlane, // copy; changed inside
                std::vector<Vertex>  &reflectedVertices,
                std::vector<Element> &reflectedElements,
                std::vector<size_t>   &vertexOrigin,
                std::vector<Isometry> &vertexIsometry)
{
    assert(onReflectionPlane.size() == 3);
    assert(onReflectionPlane[0].size() == vertices.size() &&
           onReflectionPlane[1].size() == vertices.size() &&
           onReflectionPlane[2].size() == vertices.size());

    reflectedVertices = vertices;
    reflectedElements = elements;

    // We start with the original vertices, so origins/isometries are all the
    // identity.
    vertexOrigin.assign(vertices.size(), 0);
    for (size_t i = 0; i < vertexOrigin.size(); ++i) vertexOrigin[i] = i;
    vertexIsometry.assign(vertices.size(), Isometry());

    for (size_t d = 0; d < Dim; ++d) {
        auto refl = Isometry::reflection(static_cast<Symmetry::Axis>(d));
        // We need a mapping from vertex indices of the new reflected geometry
        // we're about to create to global vertex indices.
        // All vertices except those on the reflection pane are copied.
        std::vector<size_t> globalVertexIndex(reflectedVertices.size());
        size_t numVertices = reflectedVertices.size();
        for (size_t vi = 0; vi < numVertices; ++vi) {
            if (onReflectionPlane[d].at(vi))
                globalVertexIndex[vi] = vi;
            else {
                auto v = reflectedVertices[vi];
                globalVertexIndex[vi] = reflectedVertices.size();
                v[d] *= -1;
                reflectedVertices.push_back(v);

                // Link reflected vertex back to its original in the base cell.
                vertexOrigin.push_back(vertexOrigin.at(vi));
                vertexIsometry.push_back(vertexIsometry.at(vi).compose(refl));

                // Update the onReflectionPlane info for future reflections
                for (size_t d2 = d + 1; d2 < 3; ++d2) {
                    onReflectionPlane[d2].push_back(onReflectionPlane[d2].at(vi));
                }
            }
        }
        size_t numElements = reflectedElements.size();
        for (size_t ei = 0; ei < numElements; ++ei) {
            auto re = reflectedElements[ei];
            // Reindex corner indices.
            // Note: reflection inverts the elements, so we must also permute
            // the corner indices to get positive orientation.
            // This actually matters! The inverted reflected elements cause a
            // cancellation during stiffness matrix assembly resulting in a
            // singular system.
            size_t tmp = re[0];                                                      
            re[0] = globalVertexIndex.at(re[1]);                                     
            re[1] = globalVertexIndex.at(tmp);                                       
            for (size_t d = 2; d < re.size(); ++d) re[d] = globalVertexIndex.at(re[d]);
            reflectedElements.push_back(re);                                               
        }
    }
}

#endif /* end of include guard: SNAPANDREFLECT_HH */

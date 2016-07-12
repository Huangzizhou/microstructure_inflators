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
#include <queue>
#include <ratio>
#include <MeshIO.hh>

#include <Concepts.hh>
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

// The 3D inflator is not guaranteed to place vertices precisely on the period
// cell faces. However, with the CGALClippedVolumeMesher, vertices at the
// intersection of the surface with the period cell will be placed exactly
// (they are extracted with marching squares on the cell faces).
// We can take advantage of this to segment the surface into components
// separated by period cell intersection vertices. Then we can detect and snap
// all vertices in the components lying approximately on the cell faces at
// once, which should be more robust than applying a threshold to each vertex
// independently.
// @param[inout] vertices   mesh vertices to snap
// @param[in]    m          tet mesh over "vertices"
// @param[in]    epsilon    threshold for the distance of a component's
//                          vertices to the cell boundary; if all distances of
//                          component vertices are <= this threshold, the
//                          entire component is snapped to the boundary.
template<class TMesh>
enable_if_not_models_concept_t<Concepts::TetMesh, TMesh, void>
smartSnap3D(std::vector<MeshIO::IOVertex> &/* vertices */, const TMesh &/* mesh */,
                 const Real /* epsilon */ = 1e-4) {
    throw std::runtime_error("smartSnap3D must be called on a tet mesh!");
}

template<class TMesh>
enable_if_models_concept_t<Concepts::TetMesh, TMesh, void>
smartSnap3D(std::vector<MeshIO::IOVertex> &vertices, const TMesh &mesh,
                 const Real epsilon = 1e-3) {
    const size_t nbv = mesh.numBoundaryVertices();
    std::vector<size_t> component(nbv, 0);
    std::vector<bool> visited(nbv, false);
    // Mark cell vertices as visited to segment adjacency graph
    for (auto bv : mesh.boundaryVertices()) {
        const auto &p = vertices.at(bv.volumeVertex().index()).point;
        for (size_t c = 0; c < 3; ++c) {
            if ((p[c] == 0) || (p[c] == 1)) {
                visited[bv.index()] = true;
                component[bv.index()] = 0;
                break;
            }
        }
    }

    // Segment boundary vertices into connected components
    size_t componentIdx = 0;
    for (size_t bvi = 0; bvi < nbv; ++bvi) {
        if (visited[bvi]) continue;
        ++componentIdx;

        std::queue<size_t> bfsQueue;
        bfsQueue.push(bvi);
        component[bvi] = componentIdx;
        visited[bvi] = true;

        while (!bfsQueue.empty()) {
            size_t u = bfsQueue.front();
            bfsQueue.pop();

            auto evi = mesh.boundaryVertex(u).halfEdge();
            auto eve = evi;
            do {
                assert(evi.tip().index() == u);
                size_t v = evi.tail().index();
                if (!visited[v]) {
                    visited[v] = true;
                    component[v] = componentIdx;
                    bfsQueue.push(v);
                }
            } while ((evi = evi.ccw()) != eve);
        }
    }
    const size_t nComponents = componentIdx + 1;

#if 0
    {
        // Output surface mesh marked with connected components for debugging
        ScalarField<Real> componentIndicator(nbv);
        for (size_t bvi = 0; bvi < nbv; ++bvi)
            componentIndicator[bvi] = component.at(bvi);

        std::vector<MeshIO::IOVertex > bVerts;
        std::vector<MeshIO::IOElement> bElems;
        bVerts.reserve(nbv);
        bElems.reserve(mesh.numBoundarySimplices());
        for (auto bv : mesh.boundaryVertices())
            bVerts.push_back(vertices.at(bv.volumeVertex().index()));
        for (auto bs : mesh.boundarySimplices()) {
            bElems.emplace_back(bs.vertex(0).index(),
                                bs.vertex(1).index(),
                                bs.vertex(2).index());
        }

        MSHFieldWriter writer("cellfaceComponents.msh", bVerts, bElems);
        writer.addField("component indicator", componentIndicator, DomainType::PER_NODE);
    }
#endif

    // Maximum over component vertices of distance to the 6 cell faces:
    //    min x, min y, min z, max x, max y, max z
    std::vector<std::array<Real, 6>> componentMaxDistances(nComponents, { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0} });
    for (size_t bvi = 0; bvi < nbv; ++bvi) {
        const MeshIO::IOVertex &v = vertices.at(mesh.boundaryVertex(bvi).volumeVertex().index());
        std::array<Real, 6> &cmd = componentMaxDistances.at(component[bvi]);
        // Compute dist from v to the cell faces
        for (size_t c = 0; c < 3; ++c) {
            cmd[    c] = std::max(cmd[    c], std::abs(v[c] - 0.0));
            cmd[3 + c] = std::max(cmd[3 + c], std::abs(v[c] - 1.0));
        }
    }

    const size_t NONE = std::numeric_limits<size_t>::max();
    std::vector<size_t> cellFaceForComponent(nComponents, NONE);
    for (size_t ci = 0; ci < nComponents; ++ci) {
        size_t count = 0;
        for (size_t f = 0; f < 6; ++f) {
            if (componentMaxDistances[ci][f] <= epsilon) {
                cellFaceForComponent[ci] = f;
                ++count;
            }
        }
        if (count > 1)
            std::cerr << "WARNING: ambiguous component" << std::endl;
    }

    // Snap component vertices to their respective boundaries.
    for (auto bv : mesh.boundaryVertices()) {
        size_t bvi = bv.index();
        size_t cf = cellFaceForComponent.at(component[bvi]);
        if (cf == NONE) continue;
        Real val = 0;
        if (cf >= 3) { val = 1.0; cf -= 3; }
        assert(cf < 3);
        auto &vtx = vertices.at(bv.volumeVertex().index());
        assert(std::abs(vtx.point[cf] - val) < epsilon);
        vtx[cf] = val;
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

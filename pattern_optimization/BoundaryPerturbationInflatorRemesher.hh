////////////////////////////////////////////////////////////////////////////////
// BoundaryPerturbationInflatorRemesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Remeshes the mesh created by a BoundaryPerturbationInflator. In
//      particular, short boundary segments are merged, and the interior is
//      remeshed using triangle.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/28/2015 17:56:39
////////////////////////////////////////////////////////////////////////////////
#ifndef BOUNDARYPERTURBATIONINFLATORREMESHER_HH
#define BOUNDARYPERTURBATIONINFLATORREMESHER_HH
#include <iostream>
#include <vector>
#include <queue>
#include <utils.hh>
#include "PatternOptimizationConfig.hh"
#include "BoundaryPerturbationInflator.hh"
#include <filters/extract_hole_boundaries.hh>
#include <Triangulate.h>
#include <Geometry.hh>

template<size_t N>
void remeshPerturbedShape(const BoundaryPerturbationInflator<N> &m,
                          Real maxVolume,
                          std::vector<MeshIO::IOVertex>  &outVertices,
                          std::vector<MeshIO::IOElement> &outElements);

template<>
void remeshPerturbedShape(const BoundaryPerturbationInflator<2> &m,
                          Real maxVolume,
                          std::vector<MeshIO::IOVertex>  &outVertices,
                          std::vector<MeshIO::IOElement> &outElements) {
    const auto &mesh = m.mesh();
    std::vector<Real> bdryLengths;
    size_t numBE = mesh.numBoundaryElements();
    bdryLengths.reserve(numBE);
    for (auto be : mesh.boundaryElements()) {
        auto p0 = be.vertex(0).volumeVertex().node()->p;
        auto p1 = be.vertex(1).volumeVertex().node()->p;
        bdryLengths.push_back((p1 - p0).norm());
    }
    std::vector<size_t> order;
    sortPermutation(bdryLengths, order);
    Real medianLen = bdryLengths.at(order.at(numBE / 2));

    // Build remeshed boundary edge mesh in triangulatePSLC's format
    std::vector<Point2D> bdryPts;
    std::vector<std::pair<size_t, size_t>> bdryEdges;
    // Allow computation of lengths from this format.
    auto edgeLen = [&](const std::pair<size_t, size_t> &e) { return (bdryPts.at(e.second) - bdryPts.at(e.first)).norm(); };

    // First, determine points in each hole. Assumes that hole boundary point
    // centroid lies within the hole, which should usually be true...)
    // TODO: better approach (based on winding number/shooting rays?)
    //      Don't tell triangle about holes; remove them ourselves
    //      note: this will allow triangle to insert Steiner points on the hole
    //      boundary
    std::vector<std::vector<size_t>> holeBoundaries;
    extract_hole_boundaries(mesh, holeBoundaries);
    std::vector<Point2D> holes;
    for (const auto &hb : holeBoundaries) {
        Point2D centroid(Point2D::Zero());
        for (size_t bei : hb)
            centroid += mesh.boundaryElement(bei).vertex(0).volumeVertex().node()->p;
        centroid /= hb.size();
        holes.push_back(centroid);
    }

    const auto &config = PatternOptimization::Config::get();
    Real minLength = medianLen * config.remeshMergeThreshold;
    Real maxLength = medianLen * config.remeshSplitThreshold;

    // Copy over existing boundary
    bdryPts.reserve(mesh.numBoundaryVertices());
    for (auto bv : mesh.boundaryVertices())
        bdryPts.push_back(bv.volumeVertex().node()->p);
    bdryEdges.reserve(mesh.numBoundaryElements());
    for (auto be : mesh.boundaryElements())
        bdryEdges.push_back({be.vertex(0).index(), be.vertex(1).index()});

    // Cleanup operations:
    //      1) Collapsing short edges
    //      2) Splitting  long edges
    // Perform operations in this order; edge collapse can create new long
    // edges but splitting shouldn't create short edges (for sane thresholds)

    // Collapsing removes vertices. "Sharp feature vertices" (with mesh angles
    // differing significantly from Pi) are preserved. In the absence of sharp
    // features, merges are performed in the following cases:
    //      a) Both vertices are on the same periodic boundaries.
    //         The endpoint vertices are collapsed to the edge midpoint.
    //      b) One vertex is on a subset of the other's periodic boundaries
    // Notice that collapsing will change the adjacent (next and previous) edge
    // lengths. We keep collapsing until all edges are above the collapse
    // threshold.
    // Let m represent the edge vertex appearing on a superset of the other
    // vertex' periodic faces (or vertex 0 if vertices have the same
    // membership). Collapse is executed by relocating m to the collapse
    // location, reindexing prev/next edges' vertex 1/0 to point to m, and
    // marking the boundary edge as collapsed.
    // At the end, collapsed boundary edges and unreferenced vertices are
    // deleted.

    // Mark sharp features: compute mesh angles incident on each vertex.
    std::vector<Real> meshAngles(mesh.numBoundaryVertices(), 0.0);
    for (auto t : mesh.elements()) {
        for (size_t vi = 0; vi < t.numVertices(); ++vi) {
            auto  v = t.vertex(vi);
            auto bv = v.boundaryVertex();
            if (!bv) continue;
            const auto &p0 = v.node()->p;
            //         p2 (he.tip)
            //         / \ (he)
            //        /   \ (he.tail)
            //  (v) p0-----p1
            // (in TriMesh, he_i is across from v_i)
            const auto &p1 = t.halfEdge(vi).tail().node()->p;
            const auto &p2 = t.halfEdge(vi).tip() .node()->p;
            meshAngles.at(bv.index()) += angle((p1 - p0).eval(),
                                               (p2 - p0).eval());
        }
    }

    std::vector<bool> feature(mesh.numBoundaryVertices(), false);
    for (size_t i = 0; i < meshAngles.size(); ++i) {
        if (std::abs(meshAngles[i] - M_PI) > config.remeshFeatureAngleThreshold)
            feature[i] = true;
    }

    std::queue<size_t> collapseQueue;
    for (size_t i = 0; i < order.size(); ++i) {
        if (bdryLengths.at(order[i]) < minLength)
            collapseQueue.push(order[i]);
    }
    std::vector<bool> collapsed(mesh.numBoundaryElements(), false);

    while (!collapseQueue.empty()) {
        size_t bei = collapseQueue.front();
        collapseQueue.pop();
        if (collapsed.at(bei)) continue;

        // Make sure length hasn't been increased above merge threshold
        if (bdryLengths.at(bei) > minLength) continue;
        size_t bvi0, bvi1;
        std::tie(bvi0, bvi1) = bdryEdges.at(bei);
        collapsed.at(bei) = true;

        Point2D p0 = bdryPts.at(bvi0), p1 = bdryPts.at(bvi1);
        // Can we merge 0 into 1, or 1 into 0?
        bool merge01 = m.pc().bdryVertexPeriodCellFacesPartialOrderLessEq(bvi0, bvi1);
        bool merge10 = m.pc().bdryVertexPeriodCellFacesPartialOrderLessEq(bvi1, bvi0);

        // Can't move from feature vertices.
        merge01 &= !feature.at(bvi0);
        merge10 &= !feature.at(bvi1);

        // Merge into midpoint, vertex 0, or vertex 1 location
        Point2D collapsePt;
        if (merge01 && merge10) collapsePt = (0.5 * (p0 + p1)).eval();
        else if (merge01)       collapsePt = p1;
        else if (merge10)       collapsePt = p0;
        else                    continue; // no merge

        // Call the merged vertex bvi0 unless we are merging 0 into 1
        size_t mergePtIdx = merge01 ? bvi1 : bvi0;

        // Move bv0 to the merge location.
        bdryPts.at(mergePtIdx) = collapsePt;

        // Update prev/next (uncollapsed) edges...
        auto be = mesh.boundaryElement(bei);
        auto bep = be.prev(), ben = be.next();
        // Skip over collapsed edges
        while (collapsed.at(bep.index())) bep = bep.prev();
        while (collapsed.at(ben.index())) ben = ben.next();

        size_t bepi = bep.index();
        size_t beni = ben.index();
        assert(bdryEdges.at(bepi).second == bvi0); // sanity checks
        assert(bdryEdges.at(beni).first  == bvi1); // sanity checks

        // Reindex prev/next edges' vertices to point to the merge vertex
        bdryEdges.at(bepi).second  = mergePtIdx;
        bdryEdges.at(beni).first   = mergePtIdx;

        // Update lengths, appending newly created short edges to the queue
        bdryLengths.at(bepi) = edgeLen(bdryEdges.at(bepi));
        bdryLengths.at(beni) = edgeLen(bdryEdges.at(beni));
        if (bdryLengths.at(bepi) < minLength) collapseQueue.push(bepi);
        if (bdryLengths.at(beni) < minLength) collapseQueue.push(beni);
    }

    // Delete collapsed edges
    std::vector<std::pair<size_t, size_t>> prunedBdryEdges;
    for (size_t i = 0; i < collapsed.size(); ++i)
        if (!collapsed[i]) prunedBdryEdges.push_back(bdryEdges[i]);
    bdryEdges.swap(prunedBdryEdges);
    
    // Delete unreferenced vertices, reindexing
    std::vector<bool> referenced(bdryPts.size(), false);
    for (const auto &be : bdryEdges) {
        referenced.at(be.first)  = true;
        referenced.at(be.second) = true;
    }
    std::vector<size_t> prunedBdryPtIdx(bdryPts.size(), bdryPts.size());
    std::vector<Point2D> prunedBdryPts;
    for (size_t i = 0; i < bdryPts.size(); ++i) {
        if (referenced[i]) {
            prunedBdryPtIdx[i] = prunedBdryPts.size();
            prunedBdryPts.push_back(bdryPts[i]);
        }
    }
    bdryPts.swap(prunedBdryPts);
    // Update edge endpoint indices
    for (auto &be : bdryEdges) {
        be.first  = prunedBdryPtIdx.at(be.first);
        be.second = prunedBdryPtIdx.at(be.second);
        assert(be.first < bdryPts.size());
        assert(be.second < bdryPts.size());
    }

    // Splitting introduces new vertices and can always be done.
    numBE = bdryEdges.size();
    std::vector<bool> split(numBE, false);
    for (size_t bei = 0; bei < numBE; ++bei) {
        const auto &be = bdryEdges[bei];
        if (edgeLen(be) > maxLength) {
            const auto &p0 = bdryPts.at(be.first);
            const auto &p1 = bdryPts.at(be.second);
            size_t newPt = bdryPts.size();
            // Add new vertex at midpoint
            bdryPts.push_back(0.5 * (p0 + p1));
            bdryEdges.push_back({be.first,  newPt});
            bdryEdges.push_back({newPt, be.second});
            split.push_back(false);
            split.push_back(false);
            split[bei] = true;
        }
    }
    // Delete split edges
    prunedBdryEdges.clear();
    for (size_t bei = 0; bei < bdryEdges.size(); ++bei) {
        if (!split[bei])
            prunedBdryEdges.push_back(bdryEdges[bei]);
    }
    bdryEdges.swap(prunedBdryEdges);

    {
        std::vector<MeshIO::IOVertex>  bdryOutVertices;
        std::vector<MeshIO::IOElement> bdryOutElements;
        for (const auto &p : bdryPts)
            bdryOutVertices.emplace_back(p);
        for (const auto &e : bdryEdges)
            bdryOutElements.emplace_back(e.first, e.second);
        MeshIO::save("bdry.msh", bdryOutVertices, bdryOutElements);
    }

    // Remesh the interior.
    // Q: quiet, Y: do not remesh the boundary
    triangulatePSLC(bdryPts, bdryEdges, holes, outVertices, outElements,
                    maxVolume, "QY");
}

template<>
void remeshPerturbedShape(const BoundaryPerturbationInflator<3> &/* m */,
                          Real /* maxVolume */,
                          std::vector<MeshIO::IOVertex>  &/* outVertices */,
                          std::vector<MeshIO::IOElement> &/* outElements */) {
    // 3D is a bit harder...
    throw std::runtime_error("remeshPerturbedShape<3> unimplemented.");
}

#endif /* end of include guard: BOUNDARYPERTURBATIONINFLATORREMESHER_HH */

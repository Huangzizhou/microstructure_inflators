////////////////////////////////////////////////////////////////////////////////
// WireMesh.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Pattern graph stucture with parameterized embedding.
//      Handles assignment of degrees of freedom to the graph structure based on
//      the requested symmetry and thickness parameter type (both configured by
//      template parameter)
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/26/2015 17:54:27
////////////////////////////////////////////////////////////////////////////////
#ifndef WIREMESH_HH
#define WIREMESH_HH

#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>

#include <MeshIO.hh>
#include "InflatorTypes.hh"
#include "Symmetry.hh"

enum class ThicknessType { Vertex, Edge };

template<ThicknessType thicknessType_ = ThicknessType::Vertex,
         class Symmetry_ = Symmetry::TriplyPeriodic<>>
class WireMesh {
public:
    typedef Symmetry_ PatternSymmetry;
    typedef Point3<double>            Point;
    typedef std::pair<size_t, size_t> Edge;
    // A vertex adjacent to the base cell: this is represented by the index of
    // its corresponding base cell vertex (first), and the isometry transforming
    // the base cell vertex into it (second).
    typedef std::pair<size_t, Isometry> AdjacentVertex;
    static constexpr ThicknessType thicknessType = thicknessType_;

    WireMesh(const std::string &wirePath) { load(wirePath); }

    // Embedded graph I/O (OBJ/MSH format)
    void load                  (const std::string &path);
    void save                  (const std::string &path) const;
    void saveBaseUnit          (const std::string &path) const;
    void saveReplicatedBaseUnit(const std::string &path) const;
    void saveInflationGraph    (const std::string &path, std::vector<double> params = std::vector<double>()) const;

    size_t numVertices    () const { return m_fullVertices.size(); }
    size_t numEdges       () const { return m_fullEdges   .size(); }
    size_t numBaseVertices() const { return m_baseVertices.size(); }
    size_t numBaseEdges   () const { return m_baseEdges   .size(); }

    // There is a thickness parameter for each base vertex or base edge
    // depending on thicknessType.
    size_t numThicknessParams() const {
        switch (thicknessType) {
            case ThicknessType::Edge:   return m_baseEdges.size();
            case ThicknessType::Vertex: return m_baseVertices.size();
            default: assert(false);
        }
    }
    // Positional parameters always live on the base vertices, and they are
    // determined by the base vertex positioner.
    size_t numPositionParams() const {
        size_t posParams = 0;
        for (const auto &bvp : m_baseVertexPositioners)
            posParams += bvp.numDoFs();
        return posParams;
    }

    size_t numParams() const { return numThicknessParams() + numPositionParams(); }

    // Determine the position parameters from the original embedded graph
    // positions.
    std::vector<double> defaultPositionParams() const {
        size_t np = numPositionParams();
        std::vector<double> positionParams(np);
        size_t paramOffset = 0;
        for (size_t i = 0; i < m_baseVertices.size(); ++i) {
            m_baseVertexPositioners[i].getDoFsForPoint(m_baseVertices[i],
                    &positionParams[paramOffset]);
            paramOffset += m_baseVertexPositioners[i].numDoFs();
        }
        assert(paramOffset == np);
        return positionParams;
    }

    // TODO: MAKE CONFIGURABLE.
    std::vector<double> defaultThicknessParams() const {
        return std::vector<double>(numThicknessParams(), 0.05);
    }

    // Position parameters come first, followed by thickness.
    bool isPositionParam(size_t p) const {
        if (p >= numParams()) throw std::runtime_error("Invalid parameter index.");
        return p < numPositionParams();
    }
    bool isThicknessParam(size_t p) const {
        return !isPositionParam(p);
    }

    // Position parameters come first, followed by thickness.
    std::vector<double> defaultParameters() const {
        std::vector<double> params;
        params.reserve(numParams());
        params = defaultPositionParams();
        auto dtp = defaultThicknessParams();
        params.insert(params.end(), dtp.begin(), dtp.end());
        return params;
    }

    // The inflation graph includes all vertices and edges in the base symmetry
    // unit plus adjacent edges and vertices. The adjacent subgraph is needed to
    // generate the correct inflated joint geometry at the interface nodes.
    // The "points" returned are after the positional parameters have been
    // applied, and the thickness parameters are decoded from "params" into the
    // "thicknesses" vector for convenience: these will correspond to either
    // entries in "points" or "edges" depending on thicknessType.
    template<typename Real>
    void inflationGraph(const std::vector<Real> &params,
                        std::vector<Point3<Real>> &points,
                        std::vector<Edge> &edges,
                        std::vector<Real> &thicknesses) const {
        if (params.size() != numParams())
            throw std::runtime_error("Invalid number of params.");
        edges.clear(), points.clear(), thicknesses.clear();
        edges.reserve(m_baseEdges.size() + m_adjacentEdges.size());
        points.reserve(m_baseVertices.size() + m_adjacentVertices.size());
        
        // Copy over the base graph, setting positions using params
        edges = m_baseEdges;
        size_t pOffset = 0;
        for (size_t i = 0; i < m_baseVertices.size(); ++i) {
            const auto &pos = m_baseVertexPositioners[i];
            points.push_back(pos.template getPosition<Real>(&params[pOffset]));
            pOffset += pos.numDoFs();
        }

        size_t adjVertexOffset = m_baseVertices.size();

        // Copy over the adjacent graph, transforming adjacent vertices
        // (Note that the params-repositioned vertices are transformed)
        for (const auto &ae : m_adjacentEdges)
            edges.push_back({ae.first, adjVertexOffset + ae.second});
        for (const auto &av : m_adjacentVertices)
            points.push_back(av.second.apply(points.at(av.first)));

        // params[pOffset...] are thickness parameters
        if (thicknessType == ThicknessType::Vertex) {
            // Decode into one thickness per inflation graph vertex
            thicknesses.reserve(m_baseVertices.size() + m_adjacentVertices.size());
            for (size_t i = 0; i < m_baseVertices.size(); ++i)
                thicknesses.push_back(params[pOffset + i]);
            for (size_t i = 0; i < m_adjacentVertices.size(); ++i)
                thicknesses.push_back(params[pOffset + m_adjacentVertices[i].first]);
        }
        if (thicknessType == ThicknessType::Edge) {
            // Decode into one thickness per inflation graph edge
            thicknesses.reserve(m_baseEdges.size() + m_adjacentEdges.size());
            for (size_t i = 0; i < m_baseEdges.size(); ++i)
                thicknesses.push_back(params[pOffset + i]);
            for (size_t i = 0; i < m_adjacentEdges.size(); ++i)
                thicknesses.push_back(params[pOffset + m_adjacentEdgeOrigin[i]]);
        }
    }

private:
    // All vertex/edges of the pattern
    std::vector<Point>  m_fullVertices;
    std::vector<Edge>   m_fullEdges;
    // The distinct vertices/edges of the pattern modulo symmetry
    std::vector<Point>  m_baseVertices;
    std::vector<Edge>   m_baseEdges;

    std::vector<decltype(PatternSymmetry::nodePositioner(Point()))> m_baseVertexPositioners;

    // Edges incident on one vertex in the base unit and one vertex outside.
    // The edge.first is the inside vertex, and edge.second is outside.
    // The edge.second field indexes into m_adjacentVertices
    std::vector<Edge>           m_adjacentEdges;
    std::vector<AdjacentVertex> m_adjacentVertices;
    // Which base edge originated each adjacent edge (via transormation)
    std::vector<size_t>         m_adjacentEdgeOrigin;

    // Find the base vertex within symmetry tolerance of p
    // (or throw an exception if none exists).
    size_t m_findBaseVertex(const Point &p) const;
};

#include "WireMesh.inl"

#endif /* end of include guard: WIREMESH_HH */

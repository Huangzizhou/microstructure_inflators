////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "Symmetry.hh"
#include "WireMesh.hh"
#include <json.hpp>
#include <stdexcept>
#include <memory>
#include <set>

class WireQuadMesh {
public:
    using PatternSymmetry = Symmetry::Null<>;
    using Point = WireMeshBase::Point; // Point3<double>;
    using Edge  = WireMeshBase::Edge; // std::pair<size_t, size_t>

    WireQuadMesh(const std::vector<MeshIO::IOVertex> &V, const std::vector<MeshIO::IOElement> &F, const nlohmann::json &params);

    size_t numParams() const { return m_numParams; }

    ThicknessType thicknessType() const { return m_thicknessType; }

    // Set active quad in the background mesh
    void setActiveQuad(int idx);

    // Retrieve inflation parameters for the active quad
    std::vector<double> params() const;

    // Retrieve Jacobian matrix mapping the reference square [-1,1]² to the active quad
    Eigen::Matrix3d jacobian() const;

    // Build the inflation graph for the whole quad mesh, stitching together adjacent nodes
    // (averaging stitched points' locations, thicknesses, and blending params).
    //
    // Note that the result needs to pre-warp the output point positions by the inverse of
    // the Jacobian mapping the reference square [0,1]² to the active quad. The reason for
    // that is so that `PatternSignedDistance` can be used to mesh the base cell separately
    // (the SDF-based mesher expects a base unit cell inside a specific bounding box, etc.).
    //
    // @param[in]  allParams               { Optimization variables for the active quad }
    // @param[out] stitchedPoints          { Graph vertices positions }
    // @param[out] stitchedEdges           { Graph edges indices }
    // @param[out] stitchedThicknesses     { Graph edge thicknesses }
    // @param[out] stitchedBlendingParams  { Graph vertices blending parameters }
    //
    void inflationGraph(const std::vector<double> &allParams,
        std::vector<Point>  &stitchedPoints,
        std::vector<Edge>   &stitchedEdges,
        std::vector<double> &stitchedThicknesses,
        std::vector<double> &stitchedBlendingParams) const;

private:
    nlohmann::json m_allParams;
    size_t m_numParams = 0;
    ThicknessType m_thicknessType;
};

#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#pragma once
////////////////////////////////////////////////////////////////////////////////
#include "Symmetry.hh"
#include "WireMesh.hh"
#include <json.hpp>
#include <stdexcept>
#include <memory>
#include <set>
////////////////////////////////////////////////////////////////////////////////

class WireQuadMesh {
public:
    using PatternSymmetry = Symmetry::Null<>;
    using WireMeshBasePtr = std::shared_ptr<WireMeshBase>;

    using Point = WireMeshBase::Point; // Point3<double>;
    using Edge  = WireMeshBase::Edge; // std::pair<size_t, size_t>

    WireQuadMesh(const std::vector<MeshIO::IOVertex> &V, const std::vector<MeshIO::IOElement> &F, const nlohmann::json &params);

    ThicknessType thicknessType() const { return m_thicknessType; }

    // Inflation parameters for the whole mesh (simple concatenation)
    std::vector<double> params() const;

    // Representative cell bounding box (region to be meshed)
    BBox<Point3D> boundingBox() const;

    // Build the inflation graph for the whole quad mesh, stitching together adjacent nodes
    // (averaging stitched points' locations, thicknesses, and blending params).
    //
    // @param[in]  allParams               { Inflation parameters for the whole mesh }
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
    Eigen::MatrixXd m_V;
    Eigen::MatrixXi m_F;

    std::vector<WireMeshBasePtr> m_allTopologies;
    std::vector<std::vector<double>> m_allParameters;
    std::vector<Eigen::Matrix2d> m_allJacobians; // ref square [-1,1]² to mapped parallelogram

    ThicknessType m_thicknessType = ThicknessType::Vertex;
};

#endif

#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#pragma once
////////////////////////////////////////////////////////////////////////////////
#include "BilinearMap.hh"
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

    struct MapToBaseUnit {
        BilinearMap func_;
        MapToBaseUnit(BilinearMap f = BilinearMap()) : func_(f) { }

        template<typename Real>
        Point3<Real> operator() (Point3<Real> p) const {
            return func_.apply(p[0], p[1]);
        }
    };

public:
    WireQuadMesh(const std::vector<MeshIO::IOVertex> &V, const std::vector<MeshIO::IOElement> &F, const nlohmann::json &params);

    ThicknessType thicknessType() const { return m_thicknessType; }

    // Set currently active quad
    void setActiveQuad(int idx);

    // Return wiremesh associated to the currently active quad
    const WireMeshBase &activeWireMesh() const { assert(m_activeQuad >= 0); return *m_allTopologies[m_activeQuad]; }

    // Inflation parameters for active quad
    std::vector<double> params() const { assert(m_activeQuad >= 0); return m_allParameters[m_activeQuad]; }

    MapToBaseUnit mapFunctor() const { return MapToBaseUnit(m_bilinearMap); }

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

    int m_activeQuad = -1;
    ThicknessType m_thicknessType = ThicknessType::Vertex;

    BilinearMap m_bilinearMap;

    std::vector<WireMeshBasePtr> m_allTopologies;
    std::vector<std::vector<double>> m_allParameters;
    std::vector<Eigen::Matrix3d> m_allJacobians; // ref square [-1,1]² to mapped parallelogram

    // Index of the first vertex of a quad's graph in the concatenated graph (before de-duplication)
    Eigen::VectorXi m_vertexOffset;
    Eigen::MatrixXd m_graphReducedVertices; // vertex position after de-duplication
    Eigen::VectorXi m_graphFullToReduced; // map full -> reduce vertex idx

    std::vector<std::vector<size_t>> m_stitchedVertices;
    std::vector<Edge> m_stitchedEdges;
};

#endif

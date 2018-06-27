################################################################################
# Source files
################################################################################

add_library(isosurface_inflator
        BoxIntersectionMesher.cc
        IsosurfaceInflator.cc
        MeshingOptions.cc
        MidplaneMesher.cc
        VCGSurfaceMesher.cc
        IGLSurfaceMesherMC.cc
        PostProcess.cc
        rasterize.cc

        AutomaticDifferentiation.hh
        BoxIntersection1DFeatures.hh
        BoxIntersectionMesher.hh
        CGALClippedVolumeMesher.hh
        ConvexHullTriangulation.hh
        FuzzySign.hh
        IGLSurfaceMesherMC.hh
        InflatorTypes.hh
        Isometries.hh
        IsosurfaceInflator.hh
        IsosurfaceInflatorConfig.hh
        IsosurfaceInflatorImpl.hh
        Joint.hh
        MesherBase.hh
        MeshingOptions.hh
        MidplaneMesher.hh
        NodePositioners.hh
        PaperVisualizationSDFunc.hh
        PatternSignedDistance.hh
        PostProcess.hh
        ShapeVelocity.hh
        SignedDistance.hh
        SignedDistanceRegion.hh
        SnapAndReflect.hh
        SphereConvexHull.hh
        SpherePoints.hh
        StitchedWireMesh.hh
        Symmetry.hh
        TesselateSpheres.hh
        TriangleClosestPoint.hh
        VCGSurfaceMesher.hh
        WireMesh.hh
        WireMesh.inl
        rasterize.hh
)

add_library(isosurface_inflator_cgal
        CGALClippedVolumeMesher.cc
        ConvexHullTriangulation.cc
)

target_include_directories(isosurface_inflator PUBLIC ..)

################################################################################
# Dependencies
################################################################################

target_link_libraries(isosurface_inflator_cgal
    PUBLIC
        MeshFEM
    PRIVATE
        CGAL::CGAL
)

target_link_libraries(isosurface_inflator
    PUBLIC
        MeshFEM
    PRIVATE
        isosurface_inflator_cgal
        micro::vcglib
        micro::libigl
)

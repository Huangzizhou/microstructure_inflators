################################################################################
# Source files
################################################################################

set(SOURCES
        AutomaticDifferentiation.hh
        BoxIntersection1DFeatures.hh
        BoxIntersectionMesher.cc
        BoxIntersectionMesher.hh
        CGALClippedVolumeMesher.hh
        ConvexHullTriangulation.hh
        FuzzySign.hh
        IGLSurfaceMesherMC.cc
        IGLSurfaceMesherMC.hh
        InflatorTypes.hh
        Isometries.hh
        IsosurfaceInflator.cc
        IsosurfaceInflator.hh
        IsosurfaceInflatorConfig.hh
        IsosurfaceInflatorImpl.hh
        Joint.hh
        MesherBase.hh
        MeshingOptions.cc
        MeshingOptions.hh
        MidplaneMesher.cc
        MidplaneMesher.hh
        NodePositioners.hh
        PaperVisualizationSDFunc.hh
        PatternSignedDistance.hh
        PostProcess.cc
        PostProcess.hh
        rasterize.cc
        rasterize.hh
        ShapeVelocityInterpolator.hh
        SignedDistance.hh
        SignedDistanceRegion.hh
        SnapAndReflect.hh
        SphereConvexHull.hh
        SpherePoints.hh
        StitchedWireMesh.hh
        Symmetry.hh
        TesselateSpheres.hh
        TriangleClosestPoint.hh
        TriplyPeriodicMinimalShell.hh
        VCGSurfaceMesher.cc
        VCGSurfaceMesher.hh
        WireMesh.hh
        WireMesh.inl
)

add_library(isosurface_inflator ${SOURCES})
micro_copy_headers(isosurface_inflator ${SOURCES})

add_library(isosurface_inflator_cgal
        CGALClippedVolumeMesher.cc
        ConvexHullTriangulation.cc
)

################################################################################
# Dependencies
################################################################################

target_link_libraries(isosurface_inflator_cgal
    PUBLIC
        MeshFEM
        micro::tbb
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

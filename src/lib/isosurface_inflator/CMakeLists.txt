################################################################################
# Source files
################################################################################

add_library(isosurface_inflator
        BoxIntersectionMesher.cc
        CGALClippedVolumeMesher.cc
        ConvexHullTriangulation.cc
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

target_include_directories(isosurface_inflator PUBLIC ..)

################################################################################
# Required libraries
################################################################################

# MeshFEM library
target_link_libraries(isosurface_inflator MeshFEM)

# CGAL library
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Boost 1.48 REQUIRED thread system)
target_link_libraries(isosurface_inflator CGAL::CGAL CGAL::CGAL_Core ${Boost_LIBRARIES})
target_include_directories(isosurface_inflator SYSTEM PRIVATE
    ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})

################################################################################
# Optional libraries
################################################################################

# VCG library
find_package(VCGlib)
if(VCGLIB_FOUND)
    target_link_libraries(isosurface_inflator VCGlib::core)
    target_compile_definitions(isosurface_inflator PUBLIC -DHAS_VCG)
endif()

# libigl library
find_package(LIBIGL QUIET)
if(LIBIGL_FOUND)
    target_link_libraries(isosurface_inflator igl::core)
    target_compile_definitions(isosurface_inflator PUBLIC -DHAS_LIBIGL)
endif()
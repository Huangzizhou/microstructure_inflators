################################################################################
# Common sources
################################################################################

add_library(isosurface_inflator
		BoxIntersectionMesher.cc
		CGALClippedVolumeMesher.cc
		ConvexHullTriangulation.cc
		IsosurfaceInflator.cc
		MeshingOptions.cc
		MidplaneMesher.cc
		VCGSurfaceMesher.cc)

target_include_directories(isosurface_inflator PUBLIC .)

################################################################################
# Required libraries
################################################################################

# Use C++11
set_target_properties(isosurface_inflator PROPERTIES CXX_STANDARD 11)
set_target_properties(isosurface_inflator PROPERTIES CXX_STANDARD_REQUIRED ON)

# MeshFEM library
find_package(MeshFEM REQUIRED)
target_link_libraries(isosurface_inflator MeshFEM)

# VCG library
find_package(VCGlib REQUIRED)
target_link_libraries(isosurface_inflator VCGlib::core)

# CGAL library
find_package(CGAL REQUIRED COMPONENTS Core)
set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS ON CACHE BOOL "" FORCE) # Don't mess up our build plz
include(${CGAL_USE_FILE}) # <-- Looking forwards to have imported targets with CGAL

################################################################################
# Independent binaries
################################################################################

add_executable(isosurface_cli isosurface_cli.cc)
target_link_libraries(isosurface_cli isosurface_inflator)
set_target_properties(isosurface_cli PROPERTIES CXX_STANDARD 11)

add_executable(replicate replicate.cc)
target_link_libraries(replicate isosurface_inflator)
set_target_properties(replicate PROPERTIES CXX_STANDARD 11)

add_executable(VisualizeQuantities VisualizeQuantities.cc)
target_link_libraries(VisualizeQuantities isosurface_inflator)
set_target_properties(VisualizeQuantities PROPERTIES CXX_STANDARD 11)

add_executable(PrintabilityCheck PrintabilityCheck.cc)
target_link_libraries(PrintabilityCheck isosurface_inflator)
set_target_properties(PrintabilityCheck PROPERTIES CXX_STANDARD 11)

add_executable(PrintabilityConstraints PrintabilityConstraints.cc)
target_link_libraries(PrintabilityConstraints isosurface_inflator)
set_target_properties(PrintabilityConstraints PROPERTIES CXX_STANDARD 11)

add_executable(OffsetBoundsToTranslationBounds
		OffsetBoundsToTranslationBounds.cc
		../pattern_optimization/PatternOptimizationJob.cc)
target_link_libraries(OffsetBoundsToTranslationBounds isosurface_inflator)
set_target_properties(OffsetBoundsToTranslationBounds PROPERTIES CXX_STANDARD 11)

add_executable(SignedDistance_cli SignedDistance_cli.cc)
target_link_libraries(SignedDistance_cli isosurface_inflator)
set_target_properties(SignedDistance_cli PROPERTIES CXX_STANDARD 11)

add_executable(TestWireMesh tests/TestWireMesh.cc)
target_link_libraries(TestWireMesh isosurface_inflator)
set_target_properties(TestWireMesh PROPERTIES CXX_STANDARD 11)

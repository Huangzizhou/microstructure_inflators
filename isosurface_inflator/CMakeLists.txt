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
		VCGSurfaceMesher.cc
		IGLSurfaceMesherMC.cc)

target_include_directories(isosurface_inflator PUBLIC .)

################################################################################
# Required libraries
################################################################################

# Use C++11
set_target_properties(isosurface_inflator PROPERTIES CXX_STANDARD 11)
set_target_properties(isosurface_inflator PROPERTIES CXX_STANDARD_REQUIRED ON)

# # Hide internal symbols
# set_target_properties(isosurface_inflator PROPERTIES CXX_VISIBILITY_PRESET hidden)

# MeshFEM library
find_package(MeshFEM REQUIRED)
target_link_libraries(isosurface_inflator MeshFEM)

# CGAL library
find_package(CGAL REQUIRED COMPONENTS Core)
target_link_libraries(isosurface_inflator CGAL::CGAL)

################################################################################
# Optional libraries
################################################################################

# VCG library
find_package(VCGlib)
if(VCGLIB_FOUND)
	target_link_libraries(isosurface_inflator VCGlib::core)
	target_compile_definitions(isosurface_inflator PUBLIC -DHAS_VCG)
endif()

find_package(LIBIGL QUIET)
if(LIBIGL_FOUND)
	target_compile_definitions(isosurface_inflator PUBLIC -DHAS_LIBIGL)
	target_include_directories(isosurface_inflator SYSTEM PRIVATE ${LIBIGL_INCLUDE_DIRS})
endif()

################################################################################
# Independent binaries
################################################################################

add_executable(isosurface_cli isosurface_cli.cc)
target_link_libraries(isosurface_cli isosurface_inflator)
set_target_properties(isosurface_cli PROPERTIES CXX_STANDARD 11)

add_executable(replicate replicate.cc)
target_link_libraries(replicate isosurface_inflator)
set_target_properties(replicate PROPERTIES CXX_STANDARD 11)

add_executable(tile tile.cc)
target_link_libraries(tile isosurface_inflator)
set_target_properties(tile PROPERTIES CXX_STANDARD 11)

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

add_executable(TestStitchedWireMesh tests/TestStitchedWireMesh.cc)
target_link_libraries(TestStitchedWireMesh isosurface_inflator)
set_target_properties(TestStitchedWireMesh PROPERTIES CXX_STANDARD 11)

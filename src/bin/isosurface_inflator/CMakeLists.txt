################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

add_executable(isosurface_cli isosurface_cli.cc)
target_link_libraries(isosurface_cli isosurface_inflator)

add_executable(replicate replicate.cc)
target_link_libraries(replicate isosurface_inflator)

add_executable(tile tile.cc)
target_link_libraries(tile isosurface_inflator)

add_executable(stitch_cells_cli stitch_cells_cli.cc)
target_link_libraries(stitch_cells_cli isosurface_inflator CLI11::CLI11 json)
set_target_properties(stitch_cells_cli PROPERTIES CXX_STANDARD 14)

add_executable(VisualizeQuantities VisualizeQuantities.cc)
target_link_libraries(VisualizeQuantities isosurface_inflator)

add_executable(VisualizeQuantities2D VisualizeQuantities2D.cc)
target_link_libraries(VisualizeQuantities2D isosurface_inflator)

add_executable(PrintabilityCheck PrintabilityCheck.cc)
target_link_libraries(PrintabilityCheck isosurface_inflator)

add_executable(PrintabilityConstraints PrintabilityConstraints.cc)
target_link_libraries(PrintabilityConstraints isosurface_inflator)

add_executable(SignedDistance_cli SignedDistance_cli.cc)
target_link_libraries(SignedDistance_cli isosurface_inflator)

add_executable(TestWireMesh TestWireMesh.cc)
target_link_libraries(TestWireMesh isosurface_inflator)

add_executable(TestStitchedWireMesh TestStitchedWireMesh.cc)
target_link_libraries(TestStitchedWireMesh isosurface_inflator)

find_package(LIBIGL)
if(LIBIGL_FOUND)
	add_executable(lookup_and_stitch lookup_and_stitch.cc)
	target_link_libraries(lookup_and_stitch isosurface_inflator igl::core)
endif()

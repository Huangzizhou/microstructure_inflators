################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(isosurface_cli isosurface_inflator meshfem::boost)
micro_single_app(replicate isosurface_inflator meshfem::boost)
micro_single_app(tile isosurface_inflator meshfem::boost)
micro_single_app(stitch_cells_cli isosurface_inflator CLI11::CLI11 json::json)
set_target_properties(stitch_cells_cli PROPERTIES CXX_STANDARD 14)
micro_single_app(VisualizeQuantities isosurface_inflator)
micro_single_app(VisualizeQuantities2D isosurface_inflator)
micro_single_app(PrintabilityCheck isosurface_inflator pattern_optimization)
micro_single_app(PrintabilityConstraints isosurface_inflator meshfem::boost)
micro_single_app(SignedDistance_cli isosurface_inflator meshfem::boost)
micro_single_app(TestWireMesh isosurface_inflator)
micro_single_app(TestStitchedWireMesh isosurface_inflator)
micro_single_app(TriplyPeriodicMinimalShell isosurface_inflator)
if(TARGET igl::core)
	micro_single_app(lookup_and_stitch isosurface_inflator igl::core meshfem::boost)
endif()

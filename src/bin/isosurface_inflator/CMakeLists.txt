################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(isosurface_cli micro::isosurface_inflator micro::boost)
micro_single_app(replicate micro::isosurface_inflator micro::boost)
micro_single_app(tile micro::isosurface_inflator micro::boost)
micro_single_app(stitch_cells_cli openvdb micro::isosurface_inflator CLI11::CLI11 json::json)
set_target_properties(stitch_cells_cli PROPERTIES CXX_STANDARD 14)
micro_single_app(VisualizeQuantities micro::isosurface_inflator)
micro_single_app(VisualizeQuantities2D micro::isosurface_inflator)
micro_single_app(PrintabilityConstraints micro::isosurface_inflator micro::boost)
micro_single_app(erode_cli micro::isosurface_inflator micro::boost openvdb)
micro_single_app(SignedDistance_cli micro::isosurface_inflator micro::boost)
micro_single_app(SDF_cli micro::isosurface_inflator micro::boost openvdb)
micro_single_app(TestWireMesh micro::isosurface_inflator)
micro_single_app(TestStitchedWireMesh micro::isosurface_inflator)

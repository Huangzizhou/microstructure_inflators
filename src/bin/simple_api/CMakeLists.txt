################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################

micro_single_app(homogenize MeshFEM micro::tbb CLI11::CLI11)
micro_single_app(simulate MeshFEM micro::tbb CLI11::CLI11)
micro_single_app(generate_job isosurface_inflator pattern_optimization CLI11::CLI11)
set_target_properties(generate_job PROPERTIES CXX_STANDARD 17)
micro_single_app(inflate isosurface_inflator CLI11::CLI11)
micro_single_app(inflate_and_stitch isosurface_inflator CLI11::CLI11 json::json)
set_target_properties(inflate_and_stitch PROPERTIES CXX_STANDARD 14)

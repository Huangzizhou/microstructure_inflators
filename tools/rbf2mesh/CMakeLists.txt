################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################

micro_single_app(rbf2mesh MeshFEM micro::level_set micro::inflators CLI11::CLI11)
micro_single_app(rbf_interpolation MeshFEM micro::level_set micro::inflators CLI11::CLI11)

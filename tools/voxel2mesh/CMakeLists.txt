################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################

micro_single_app(voxel2mesh MeshFEM micro::voxel_inflator micro::inflators CLI11::CLI11)

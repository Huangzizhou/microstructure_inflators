################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

add_executable(homogenize homogenize.cc)
target_link_libraries(homogenize MeshFEM CLI11::CLI11)

add_executable(simulate simulate.cc)
target_link_libraries(simulate MeshFEM CLI11::CLI11)

add_executable(inflate inflate.cc)
target_link_libraries(inflate isosurface_inflator CLI11::CLI11)

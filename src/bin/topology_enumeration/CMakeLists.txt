################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

add_library(intersection_check intersection_check.cc)
target_link_libraries(intersection_check isosurface_inflator MeshFEM micro::cgal)

################################################################################
# Independent binaries
################################################################################

add_executable(enumerate enumerate.cc)
target_link_libraries(enumerate MeshFEM intersection_check)

add_executable(enumerate_2D enumerate.cc)
target_compile_definitions(enumerate_2D PUBLIC -DDIM=2)
target_link_libraries(enumerate_2D MeshFEM intersection_check)

add_executable(match match.cc)
target_link_libraries(match MeshFEM)

add_executable(sort sort.cc)
target_link_libraries(sort MeshFEM)

add_executable(pattern_info pattern_info.cc)
target_link_libraries(pattern_info MeshFEM isosurface_inflator)

# MeshFEM library
find_package(MeshFEM REQUIRED)

# CGAL library
find_package(CGAL REQUIRED COMPONENTS Core)

add_library(intersection_check intersection_check.cc)
target_link_libraries(intersection_check MeshFEM)
set_target_properties(intersection_check PROPERTIES CXX_STANDARD 11)
set_target_properties(intersection_check PROPERTIES CXX_STANDARD_REQUIRED ON)
target_link_libraries(intersection_check CGAL::CGAL)
target_include_directories(intersection_check SYSTEM PRIVATE ${CGAL_INCLUDE_DIRS})

################################################################################
# Independent binaries
################################################################################
add_executable(enumerate enumerate.cc)
target_link_libraries(enumerate MeshFEM intersection_check)
set_target_properties(enumerate PROPERTIES CXX_STANDARD 11)

add_executable(enumerate_2D enumerate.cc)
target_compile_definitions(enumerate_2D PUBLIC -DDIM=2)
target_link_libraries(enumerate_2D MeshFEM intersection_check)
set_target_properties(enumerate_2D PROPERTIES CXX_STANDARD 11)

add_executable(match match.cc)
target_link_libraries(match MeshFEM)
set_target_properties(match PROPERTIES CXX_STANDARD 11)

add_executable(sort sort.cc)
target_link_libraries(sort MeshFEM)
set_target_properties(sort PROPERTIES CXX_STANDARD 11)

add_executable(pattern_info pattern_info.cc)
target_link_libraries(pattern_info MeshFEM)
set_target_properties(pattern_info PROPERTIES CXX_STANDARD 11)

add_library(intersection_check intersection_check.cc)

# MeshFEM library
target_link_libraries(intersection_check MeshFEM isosurface_inflator)

# CGAL library
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Boost 1.48 REQUIRED thread system)
target_link_libraries(intersection_check CGAL::CGAL CGAL::CGAL_Core ${Boost_LIBRARIES})
target_include_directories(intersection_check SYSTEM PRIVATE
    ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})

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

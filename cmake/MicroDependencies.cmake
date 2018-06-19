# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

### Configuration
set(MICRO_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(MICRO_EXTERNAL "${MICRO_ROOT}/3rdparty")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(MicroDownloadExternal)

################################################################################
# Required libraries
################################################################################

# MeshFEM library
if(NOT TARGET MeshFEM)
    add_subdirectory(${MICRO_EXTERNAL}/MeshFEM MeshFEM)
endif()

# CLI11 library
if(NOT TARGET CLI11::CLI11)
    add_subdirectory(${MICRO_EXTERNAL}/CLI11)
endif()

# CGAL library
if(NOT TARGET micro::cgal)
    add_library(micro_cgal INTERFACE)
    find_package(CGAL REQUIRED COMPONENTS Core)
    find_package(Boost 1.48 REQUIRED COMPONENTS thread system)
    target_link_libraries(micro_cgal INTERFACE CGAL::CGAL CGAL::CGAL_Core ${Boost_LIBRARIES})
    target_include_directories(micro_cgal SYSTEM INTERFACE
        ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
    add_library(micro::cgal ALIAS micro_cgal)
endif()

# Accelerate framework
if(NOT TARGET micro::accelerate)
    if(APPLE)
        find_library(AccelerateFramework Accelerate)
        add_library(micro_accelerate INTERFACE)
        target_link_libraries(micro_accelerate INTERFACE ${AccelerateFramework})
        add_library(micro::accelerate ALIAS micro_accelerate)
    else()
        add_library(micro::accelerate INTERFACE IMPORTED)
    endif()
endif()

################################################################################
# Optional libraries
################################################################################

# Clipper library
if(NOT TARGET micro::clipper)
    find_package(Clipper QUIET)
    if(CLIPPER_FOUND)
        add_library(micro_clipper INTERFACE)
        target_include_directories(micro_clipper SYSTEM INTERFACE ${CLIPPER_INCLUDE_DIRS})
        target_link_libraries(micro_clipper INTERFACE ${CLIPPER_LIBRARIES})
        target_compile_definitions(micro_clipper INTERFACE -DHAS_CLIPPER)
        add_library(micro::clipper ALIAS micro_clipper)
    else()
        message(STATUS "Clipper not found; disabling Luigi's Inflator")
        add_library(micro::clipper INTERFACE IMPORTED)
    endif()
endif()

# Ceres library
if(NOT TARGET micro::ceres)
    find_package(Ceres QUIET)
    if(CERES_FOUND)
        add_library(micro_ceres INTERFACE)
        target_include_directories(micro_ceres SYSTEM INTERFACE ${CERES_INCLUDE_DIRS})
        target_link_libraries(micro_ceres INTERFACE ${CERES_LIBRARIES})
        target_compile_definitions(micro_ceres INTERFACE -DHAS_CERES)
        add_library(micro::ceres ALIAS micro_ceres)
    else()
        message(STATUS "Google's ceres-solver not found; levenberg-marquardt disabled")
        add_library(micro::ceres INTERFACE IMPORTED)
    endif()
endif()

# Dlib Library
if(NOT TARGET micro::dlib)
    find_package(Dlib QUIET)
    if(DLIB_FOUND)
        # Complains when Dlib is compiled in Release and the main app is compiled in Debug
        # target_include_directories(optimizers SYSTEM PUBLIC ${DLIB_INCLUDE_DIR})
        # target_link_libraries(optimizers ${DLIB_LIBRARIES})
        # target_compile_definitions(optimizers PUBLIC -DHAS_DLIB)
        add_library(micro::dlib INTERFACE IMPORTED)
    else()
        message(STATUS "DLib not found; bfgs optimizer is disabled")
        add_library(micro::dlib INTERFACE IMPORTED)
    endif()
endif()

# Knitro library
if(NOT TARGET micro::knitro)
    find_package(Knitro QUIET)
    if(KNITRO_FOUND)
        add_library(micro_knitro INTERFACE)
        target_include_directories(micro_knitro SYSTEM INTERFACE ${KNITRO_INCLUDE_DIRS})
        target_link_libraries(micro_knitro INTERFACE ${KNITRO_LIBRARIES})
        target_compile_definitions(micro_knitro INTERFACE -DHAS_KNITRO)
        add_library(micro::knitro ALIAS micro_knitro)
    else()
        message(STATUS "Knitro not found; active_set disabled")
        add_library(micro::knitro INTERFACE IMPORTED)
    endif()
endif()

# libigl library
if(NOT TARGET micro::libigl)
    find_package(LIBIGL QUIET)
    if(LIBIGL_FOUND)
        add_library(micro_libigl INTERFACE)
        target_link_libraries(micro_libigl INTERFACE igl::core)
        target_compile_definitions(micro_libigl INTERFACE -DHAS_LIBIGL)
        add_library(micro::libigl ALIAS micro_libigl)
    else()
        add_library(micro::libigl INTERFACE IMPORTED)
    endif()
endif()

# NLopt library
if(NOT TARGET micro::nlopt)
    find_package(NLopt QUIET)
    if(NLOPT_FOUND)
        add_library(micro_nlopt INTERFACE)
        target_link_libraries(micro_nlopt INTERFACE ${NLOPT_LIBRARIES})
        target_compile_definitions(micro_nlopt INTERFACE -DHAS_NLOPT)
        add_library(micro::nlopt ALIAS micro_nlopt)
    else()
        message(STATUS "NLopt not found; slsqp disabled")
        add_library(micro::nlopt INTERFACE IMPORTED)
    endif()
endif()

# PyMesh Wires library
if(NOT TARGET micro::pymesh)
    find_package(PyMesh QUIET)
    if(PYMESH_FOUND)
        add_library(micro_pymesh INTERFACE)
        target_link_libraries(micro_pymesh INTERFACE PyMesh::core PyMesh::wires)
        target_compile_definitions(micro_pymesh INTERFACE -DHAS_PYMESH)
        add_library(micro::pymesh ALIAS micro_pymesh)
    else()
        message(STATUS "PyMesh not found; disabling James' Inflator")
        add_library(micro::pymesh INTERFACE IMPORTED)
    endif()
endif()

# VCG library
if(NOT TARGET micro::vcglib)
    find_package(VCGlib QUIET)
    if(VCGLIB_FOUND)
        add_library(micro_vcglib INTERFACE)
        target_link_libraries(micro_vcglib INTERFACE VCGlib::core)
        # target_compile_definitions(inflators PUBLIC -DHAS_VCGLIB)
        add_library(micro::vcglib ALIAS micro_vcglib)
    else()
        message(STATUS "VCGLib not found; disabling Luigi's Inflator")
        add_library(micro::vcglib INTERFACE IMPORTED)
    endif()
endif()

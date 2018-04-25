################################################################################
cmake_minimum_required(VERSION 3.1)
project(Microstructures_External)
################################################################################

include(ExternalProject)

################################################################################
# Platform-specific options
################################################################################

cmake_host_system_information(RESULT MY_HOSTNAME QUERY HOSTNAME)
if(MY_HOSTNAME MATCHES "log-*")
    set(DISABLE_CXX11_ABI_DEFAULT ON)
else()
    set(DISABLE_CXX11_ABI_DEFAULT OFF)
endif()
option(DISABLE_CXX11_ABI_GCC "Don't use GCC's new C++11 ABI; needed to prevent linker errors with libraries compiled with the old ABI" ${DISABLE_CXX11_ABI_DEFAULT})

if(DISABLE_CXX11_ABI_GCC)
  set(CMAKE_EXTRA_ARGS "-DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0")
  set(MAKE_EXTRA_ARGS "_GLIBCXX_USE_CXX11_ABI=0")
else()
  set(CMAKE_EXTRA_ARGS "")
  set(MAKE_EXTRA_ARGS "")
endif()

################################################################################
# Ceres Solver
################################################################################

set(CERES_VERSION "1.13.0")
set(CERES_URL http://ceres-solver.org/ceres-solver-${CERES_VERSION}.tar.gz)
set(CERES_MD5 cd568707571c92af3d69c1eb28d63d72)
set(CERES_DIR "${CMAKE_CURRENT_LIST_DIR}/ceres")

ExternalProject_Add(ceres_external
    PREFIX "${CERES_DIR}"
    URL ${CERES_URL}
    URL_MD5 ${CERES_MD5}
    LIST_SEPARATOR |
    CMAKE_ARGS
        "-DMINIGLOG:BOOL=OFF"
        "-DCMAKE_INSTALL_PREFIX=${CERES_DIR}/"
        "-DEIGEN_INCLUDE_DIR_HINTS:PATH=$ENV{EIGEN_INC}"
        "-DSUITESPARSE_INCLUDE_DIR_HINTS:LIST=$ENV{SUITESPARSE_INC}|$ENV{TBB_INC}"
        "-DSUITESPARSE_LIBRARY_DIR_HINTS:LIST=$ENV{SUITESPARSE_LIB}|$ENV{TBB_LIB}"
        "-DGLOG_INCLUDE_DIR_HINTS:PATH=$ENV{GLOG_INC}"
        "-DGLOG_LIBRARY_DIR_HINTS:PATH=$ENV{GLOG_LIB}"
        "${CMAKE_EXTRA_ARGS}"
)

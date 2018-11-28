################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(MICRO_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(MICRO_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(micro_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${MICRO_EXTERNAL}/${name}
        DOWNLOAD_DIR ${MICRO_EXTERNAL}/.cache/${name}
        QUIET
        ${MICRO_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## MeshFEM
function(micro_download_meshfem)
    micro_download_project(MeshFEM
        GIT_REPOSITORY git@github.com:geometryprocessing/MeshFEM.git
        GIT_TAG        a6673722148332c8a9b1799432f81c02624cc367
    )
endfunction()

## TBB
function(micro_download_tbb)
    micro_download_project(tbb
        GIT_REPOSITORY https://github.com/01org/tbb.git
        GIT_TAG        2019_U1
        # GIT_TAG        08b4341a1893a72656467e96137f1f99d0112547
    )
endfunction()

## CGAL
function(micro_download_cgal)
    micro_download_project(cgal
        URL     https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.12/CGAL-4.12.tar.xz
        URL_MD5 b12fd24dedfa889a04abfaea565a88bd
    )
endfunction()

## Catch2
function(micro_download_catch)
    micro_download_project(Catch2
        URL     https://github.com/catchorg/Catch2/archive/v2.3.0.tar.gz
        URL_MD5 1fc90ff3b7b407b83057537f4136489e
    )
endfunction()

## CLI11
function(micro_download_cli11)
    micro_download_project(CLI11
        URL     https://github.com/CLIUtils/CLI11/archive/v1.6.0.tar.gz
        URL_MD5 c8e3dc70e3b7ebf6b01f618f7cdcc85f
    )
endfunction()

## nlopt
function(micro_download_nlopt)
    micro_download_project(nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG        8da50d516fe0a0da42c06b66329b02cf6e44236e
    )
endfunction()

## libigl
function(micro_download_libigl)
    micro_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        75d60e40a8edc6868571fbdca2e74f97d5dddab8
    )
endfunction()

## Sanitizers
function(micro_download_sanitizers)
    micro_download_project(sanitizers-cmake
        GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
        GIT_TAG        6947cff3a9c9305eb9c16135dd81da3feb4bf87f
    )
endfunction()

## Cotire
function(micro_download_cotire)
    micro_download_project(cotire
        GIT_REPOSITORY https://github.com/sakra/cotire.git
        GIT_TAG        391bf6b7609e14f5976bd5247b68d63cbf8d4d12
    )
endfunction()

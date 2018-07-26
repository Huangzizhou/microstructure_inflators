################################################################################
include(DownloadProject)

# Shortcut function
function(micro_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${MICRO_EXTERNAL}/${name}
        DOWNLOAD_DIR ${MICRO_EXTERNAL}/.cache/${name}
        ${ARGN}
    )
endfunction()

################################################################################

## TBB
function(micro_download_tbb)
    micro_download_project(tbb
        GIT_REPOSITORY https://github.com/wjakob/tbb.git
        GIT_TAG        4c3ffe5a5f37addef0dd6283c74c4402a3b4ebc9
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
        # GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        # GIT_TAG        58444833fdafdfc2046f8fdfdc5a40bdb486b017
        GIT_REPOSITORY https://github.com/jdumas/nlopt.git
        GIT_TAG        269c997f2703915c0d32c2534e6734a51b51514e
    )
endfunction()

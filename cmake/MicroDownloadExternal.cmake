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

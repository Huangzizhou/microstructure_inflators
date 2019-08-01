################################################################################
# Source files
################################################################################

set(SOURCES
        RBF.cc
        RBF.hh
)

micro_add_library(level_set ${SOURCES})
micro_copy_headers(level_set ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_level_set
    PUBLIC
        MeshFEM
        meshfem::boost
    PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(micro_level_set)
endif()
################################################################################
# Source files
################################################################################

set(SOURCES
        VoxelSDF.cc
        VoxelSDF.hh
)

micro_add_library(voxel_inflator ${SOURCES})
micro_copy_headers(voxel_inflator ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_voxel_inflator
    PUBLIC
        MeshFEM
        meshfem::boost
        micro::hex_inflator
        micro::isosurface_inflator
    PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(micro_voxel_inflator)
endif()
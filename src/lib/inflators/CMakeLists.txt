################################################################################
# Source files
################################################################################

set(SOURCES
        MakeInflator.cc
        Inflator.hh
        MakeInflator.hh
        wrappers/BoundaryPerturbationInflator.cc
        wrappers/BoundaryPerturbationInflator.hh
        wrappers/ConstrainedIsoinflator.cc
        wrappers/ConstrainedIsoinflator.hh
        wrappers/EqualityConstrainedInflator.cc
        wrappers/EqualityConstrainedInflator.hh
        wrappers/HexaPillarsInflator.cc
        wrappers/HexaPillarsInflator.hh
        wrappers/IsoinflatorWrapper.cc
        wrappers/IsoinflatorWrapper.hh
        wrappers/JamesInflatorWrapper.cc
        wrappers/JamesInflatorWrapper.hh
        wrappers/LuigiInflatorWrapper.cc
        wrappers/LuigiInflatorWrapper.hh
        wrappers/VoxelsInflator.cc
        wrappers/VoxelsInflator.hh
        wrappers/RBFInflator.cc
        wrappers/RBFInflator.hh
        wrappers/ConstrainedInflator.hh
)

micro_add_library(inflators ${SOURCES})
micro_copy_headers(inflators ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_inflators
    PUBLIC
        MeshFEM
        meshfem::boost
        micro::hex_inflator
        micro::isosurface_inflator
        micro::voxel_inflator
        micro::pymesh
        micro::vcglib
        micro::clipper
    PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(micro_inflators)
endif()

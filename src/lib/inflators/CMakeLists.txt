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
)

add_library(inflators ${SOURCES})
micro_copy_headers(inflators ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(inflators
    PUBLIC
        MeshFEM
        hex_inflator
        isosurface_inflator
        micro::pymesh
        micro::vcglib
        micro::clipper
    PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(inflators)
endif()

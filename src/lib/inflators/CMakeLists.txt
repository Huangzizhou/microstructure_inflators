################################################################################
# Source files
################################################################################

add_library(inflators
        MakeInflator.cc
        wrappers/EqualityConstrainedInflator.cc
        wrappers/IsoinflatorWrapper.cc
        wrappers/JamesInflatorWrapper.cc
        wrappers/LuigiInflatorWrapper.cc
        wrappers/BoundaryPerturbationInflator.cc
        wrappers/ConstrainedIsoinflator.cc
        wrappers/HexaPillarsInflator.cc
)

target_include_directories(inflators PUBLIC ..)

################################################################################
# Dependencies
################################################################################

target_link_libraries(inflators
        MeshFEM
        micro::pymesh
        micro::vcglib
        micro::clipper
)

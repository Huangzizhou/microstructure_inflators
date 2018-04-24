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
# Required libraries
################################################################################

# MeshFEM library
target_link_libraries(inflators MeshFEM json)

################################################################################
# Optional libraries
################################################################################

# PyMesh Wires library
find_package(PyMesh QUIET)
if(PYMESH_FOUND)
    target_link_libraries(inflators PyMesh::core PyMesh::wires)
    target_compile_definitions(inflators PUBLIC -DHAS_PYMESH)
else()
    message(STATUS "PyMesh not found; disabling James' Inflator")
endif()

# VCG library
find_package(VCGlib QUIET)
if(VCGLIB_FOUND)
    target_link_libraries(inflators VCGlib::core)
    # target_compile_definitions(inflators PUBLIC -DHAS_VCGLIB)
else()
    message(STATUS "VCGLib not found; disabling Luigi's Inflator")
endif()

# Clipper library
find_package(Clipper QUIET)
if(CLIPPER_FOUND)
    target_include_directories(inflators SYSTEM PUBLIC ${CLIPPER_INCLUDE_DIRS})
    target_link_libraries(inflators ${CLIPPER_LIBRARIES})
    target_compile_definitions(inflators PUBLIC -DHAS_CLIPPER)
else()
    message(STATUS "Clipper not found; disabling Luigi's Inflator")
endif()

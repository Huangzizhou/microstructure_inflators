################################################################################
# Common sources
################################################################################

add_library(pattern_optimization PatternOptimizationJob.cc)
target_include_directories(pattern_optimization PUBLIC .)
set_target_properties(pattern_optimization PROPERTIES CXX_STANDARD 11)

# isosurface_inflator
target_link_libraries(pattern_optimization isosurface_inflator)

################################################################################

add_library(inflators
        MakeInflator.cc
        inflators/EqualityConstrainedInflator.cc
        inflators/IsoinflatorWrapper.cc
        inflators/JamesInflatorWrapper.cc
        inflators/LuigiInflatorWrapper.cc
		inflators/BoundaryPerturbationInflator.cc
		inflators/HexaPillarsInflator.cpp)
set_target_properties(inflators PROPERTIES CXX_STANDARD 11)

# MeshFEM library
target_link_libraries(inflators MeshFEM json)

################################################################################
# Optional libraries for inflation
################################################################################

# PyMesh Wires library
#find_package(PyMesh QUIET)
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
    target_compile_definitions(inflators PUBLIC -DHAS_VCGLIB)
else()
    message(STATUS "VCGLib not found; disabling Luigi's Inflator")
endif()

# Clipper library
#find_package(Clipper QUIET)
if(CLIPPER_FOUND)
    target_include_directories(inflators SYSTEM PUBLIC ${CLIPPER_INCLUDE_DIRS})
    target_link_libraries(     inflators        ${CLIPPER_LIBRARIES})
    target_compile_definitions(inflators PUBLIC -DHAS_CLIPPER)
else()
    message(STATUS "Clipper not found; disabling Luigi's Inflator")
endif()

################################################################################
add_library(optimizers
        optimizers/ceres.cc
        optimizers/dlib.cc
        optimizers/gradient_descent.cc
        optimizers/nlopt.cc
        optimizers/knitro.cc)
set_target_properties(optimizers PROPERTIES CXX_STANDARD 11)

# MeshFEM library
target_link_libraries(optimizers MeshFEM)

# Accelerate framework
if (APPLE)
    find_library(AccelerateFramework Accelerate)
    target_link_libraries(optimizers ${AccelerateFramework})
endif (APPLE)

################################################################################
# Optional optimization libraries
################################################################################

# Dlib Library
find_package(Dlib QUIET)
if(DLIB_FOUND)
    target_include_directories(optimizers SYSTEM PUBLIC ${DLIB_INCLUDE_DIR})
    target_link_libraries(optimizers ${DLIB_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_DLIB)
else()
    message(STATUS "DLib not found; bfgs optimizer is disabled")
endif()

# Ceres library
find_package(Ceres QUIET)
if(CERES_FOUND)
    target_include_directories(optimizers SYSTEM PUBLIC ${CERES_INCLUDE_DIRS})
    target_link_libraries(optimizers ${CERES_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_CERES)
else()
    message(STATUS "Google's ceres-solver not found; levenberg-marquardt disabled")
endif()

# NLopt library
find_package(NLopt QUIET)
if(NLOPT_FOUND)
    target_include_directories(optimizers SYSTEM PUBLIC ${NLOPT_INCLUDE_DIR})
    target_link_libraries(optimizers ${NLOPT_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_NLOPT)
else()
    message(STATUS "NLopt not found; slsqp disabled")
endif()

# Knitro library
find_package(Knitro QUIET)
if(KNITRO_FOUND)
    target_include_directories(optimizers SYSTEM PUBLIC ${KNITRO_INCLUDE_DIRS})
    target_link_libraries(optimizers ${KNITRO_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_KNITRO)
else()
    message(STATUS "Knitro not found; active_set disabled")
endif()

################################################################################
# Independent binaries
################################################################################

add_executable(PatternOptimization_cli PatternOptimization_cli.cc)
target_link_libraries(PatternOptimization_cli pattern_optimization optimizers inflators)
set_target_properties(PatternOptimization_cli PROPERTIES CXX_STANDARD 11)

add_executable(GenIsosurfaceJob GenIsosurfaceJob.cc)
target_link_libraries(GenIsosurfaceJob pattern_optimization)
set_target_properties(GenIsosurfaceJob PROPERTIES CXX_STANDARD 14)

# add_executable(RandomJob RandomJob.cc)
# target_link_libraries(RandomJob pattern_optimization)
# set_target_properties(RandomJob PROPERTIES CXX_STANDARD 11)

# add_executable(SpacedJob SpacedJob.cc)
# target_link_libraries(SpacedJob pattern_optimization)
# set_target_properties(SpacedJob PROPERTIES CXX_STANDARD 11)

# add_executable(CreateJob CreateJob.cc)
# target_link_libraries(CreateJob pattern_optimization)
# set_target_properties(CreateJob PROPERTIES CXX_STANDARD 11)

# add_executable(GradientComponentValidation GradientComponentValidation.cc)
# target_link_libraries(GradientComponentValidation pattern_optimization)
# set_target_properties(GradientComponentValidation PROPERTIES CXX_STANDARD 11)

# add_executable(Inflator_cli Inflator_cli.cc)
# target_link_libraries(Inflator_cli pattern_optimization)
# set_target_properties(Inflator_cli PROPERTIES CXX_STANDARD 11)

# add_executable(Inflator2D_cli Inflator_cli.cc)
# target_link_libraries(Inflator2D_cli pattern_optimization)
# target_compile_definitions(Inflator2D_cli PUBLIC -DDIMENSIONS=2)
# set_target_properties(Inflator2D_cli PROPERTIES CXX_STANDARD 11)

# add_executable(ParamRegion_cli ParamRegion_cli.cc)
# target_link_libraries(ParamRegion_cli pattern_optimization)
# set_target_properties(ParamRegion_cli PROPERTIES CXX_STANDARD 11)

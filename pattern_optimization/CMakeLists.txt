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
		inflators/BoundaryPerturbationInflator.cc)
set_target_properties(inflators PROPERTIES CXX_STANDARD 11)

# PyMesh Wires library
find_package(PyMesh REQUIRED)
target_link_libraries(inflators PyMesh::core PyMesh::wires)

# MeshFEM library
target_link_libraries(inflators MeshFEM)

# VCG library
find_package(VCGlib REQUIRED)
target_link_libraries(inflators VCGlib::core)

################################################################################

add_library(optimizers
		optimizers/ceres.cc
		optimizers/dlib.cc
		optimizers/gradient_descent.cc
		optimizers/nlopt.cc)
set_target_properties(optimizers PROPERTIES CXX_STANDARD 11)

# MeshFEM library
target_link_libraries(optimizers MeshFEM)

# NLopt library
find_package(NLopt REQUIRED)
target_include_directories(optimizers PUBLIC ${NLOPT_INCLUDE_DIR})
target_link_libraries(optimizers ${NLOPT_LIBRARIES})

################################################################################

# CGAL library
find_package(CGAL REQUIRED COMPONENTS Core)
set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS ON CACHE BOOL "" FORCE) # Don't mess up our build plz
include(${CGAL_USE_FILE}) # <-- Looking forwards to have imported targets with CGAL

################################################################################
# Independent binaries
################################################################################

add_executable(PatternOptimization_cli PatternOptimization_cli.cc)
target_link_libraries(PatternOptimization_cli pattern_optimization optimizers inflators)
set_target_properties(PatternOptimization_cli PROPERTIES CXX_STANDARD 11)

add_executable(GenIsosurfaceJob GenIsosurfaceJob.cc)
target_link_libraries(GenIsosurfaceJob pattern_optimization)
set_target_properties(GenIsosurfaceJob PROPERTIES CXX_STANDARD 11)

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

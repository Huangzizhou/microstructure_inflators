################################################################################
# Common sources
################################################################################

# None

################################################################################
# Independent binaries
################################################################################

# add_executable(DiscreteShapeDerivativeValidation_cli DiscreteShapeDerivativeValidation_cli.cc)
# target_link_libraries(DiscreteShapeDerivativeValidation_cli MeshFEM)
# set_target_properties(DiscreteShapeDerivativeValidation_cli PROPERTIES CXX_STANDARD 11)

add_executable(MacroMicroStress_cli MacroMicroStress_cli.cc)
target_link_libraries(MacroMicroStress_cli MeshFEM)
set_target_properties(MacroMicroStress_cli PROPERTIES CXX_STANDARD 11)

add_executable(WCS_cli WCS_cli.cc)
target_link_libraries(WCS_cli MeshFEM)
set_target_properties(WCS_cli PROPERTIES CXX_STANDARD 11)

add_executable(WCSOptimization_cli WCSOptimization_cli.cc)
target_link_libraries(WCSOptimization_cli pattern_optimization optimizers inflators json)
set_target_properties(WCSOptimization_cli PROPERTIES CXX_STANDARD 11)

# add_executable(WCSRemeshingOptimization_cli WCSRemeshingOptimization_cli.cc)
# target_link_libraries(WCSRemeshingOptimization_cli pattern_optimization)
# set_target_properties(WCSRemeshingOptimization_cli PROPERTIES CXX_STANDARD 11)

# making CLION happy
#FILE (GLOB_RECURSE clion_all_headers ${CMAKE_SOURCE_DIR}/*.hh)
#FILE (GLOB_RECURSE clion_all_sources ${CMAKE_SOURCE_DIR}/*.cc)
#add_executable(makingClionHappy EXCLUDE_FROM_ALL ${clion_all_headers} ${clion_all_sources})

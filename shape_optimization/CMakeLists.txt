################################################################################
# Common sources
################################################################################

# None

################################################################################
# Independent binaries
################################################################################

add_executable(ShapeOptimization_cli ShapeOptimization_cli.cc)
target_link_libraries(ShapeOptimization_cli pattern_optimization optimizers inflators json)
set_target_properties(ShapeOptimization_cli PROPERTIES CXX_STANDARD 11)

add_executable(ShapeDerivativeValidation_cli ShapeDerivativeValidation_cli.cc)
target_link_libraries(ShapeDerivativeValidation_cli pattern_optimization inflators MeshFEM)
set_target_properties(ShapeDerivativeValidation_cli PROPERTIES CXX_STANDARD 11)

add_executable(GenerateParametersMask GenerateParametersMask.cc)
target_link_libraries(GenerateParametersMask pattern_optimization inflators json)
set_target_properties(GenerateParametersMask PROPERTIES CXX_STANDARD 11)
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

add_executable(GenerateParametersMask_cli GenerateParametersMask_cli.cc)
target_link_libraries(GenerateParametersMask_cli pattern_optimization inflators)
set_target_properties(GenerateParametersMask_cli PROPERTIES CXX_STANDARD 11)

add_executable(IsosurfaceShapeDerivativeValidation_cli IsosurfaceShapeDerivativeValidation_cli.cc)
target_link_libraries(IsosurfaceShapeDerivativeValidation_cli pattern_optimization inflators MeshFEM)
set_target_properties(IsosurfaceShapeDerivativeValidation_cli PROPERTIES CXX_STANDARD 11)

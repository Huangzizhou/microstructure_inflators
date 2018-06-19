################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(ShapeOptimization_cli pattern_optimization optimizers inflators json::json)
micro_single_app(ShapeDerivativeValidation_cli pattern_optimization inflators MeshFEM)
micro_single_app(GenerateParametersMask_cli pattern_optimization inflators)
micro_single_app(IsosurfaceShapeDerivativeValidation_cli pattern_optimization inflators MeshFEM)
micro_single_app(VolumeShapeDerivativeValidation_cli pattern_optimization inflators MeshFEM)

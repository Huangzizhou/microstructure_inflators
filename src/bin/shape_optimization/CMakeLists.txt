################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(ShapeOptimization_cli micro::pattern_optimization micro::optimizers micro::inflators json::json meshfem::boost)
micro_single_app(ShapeDerivativeValidation_cli micro::pattern_optimization micro::inflators MeshFEM meshfem::boost)
micro_single_app(GenerateParametersMask_cli micro::pattern_optimization micro::inflators meshfem::boost)
micro_single_app(IsosurfaceShapeDerivativeValidation_cli micro::pattern_optimization micro::inflators MeshFEM meshfem::boost)
micro_single_app(VolumeShapeDerivativeValidation_cli micro::pattern_optimization micro::inflators MeshFEM meshfem::boost)

################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(DiscreteShapeDerivativeValidation_cli MeshFEM pattern_optimization inflators micro::tbb)
micro_single_app(MacroMicroStress_cli pattern_optimization MeshFEM micro::tbb)
micro_single_app(WCS_cli pattern_optimization MeshFEM micro::tbb)
micro_single_app(WCSOptimization_cli pattern_optimization optimizers inflators json::json micro::tbb)
# micro_single_app(WCSRemeshingOptimization_cli pattern_optimization)

################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(PatternOptimization_cli pattern_optimization optimizers inflators)
micro_single_app(GenIsosurfaceJob pattern_optimization)
set_target_properties(GenIsosurfaceJob PROPERTIES CXX_STANDARD 14)
micro_single_app(OffsetBoundsToTranslationBounds isosurface_inflator pattern_optimization)

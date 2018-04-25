################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

add_executable(PatternOptimization_cli PatternOptimization_cli.cc)
target_link_libraries(PatternOptimization_cli pattern_optimization optimizers inflators)

add_executable(GenIsosurfaceJob GenIsosurfaceJob.cc)
target_link_libraries(GenIsosurfaceJob pattern_optimization)
set_target_properties(GenIsosurfaceJob PROPERTIES CXX_STANDARD 14)

add_executable(OffsetBoundsToTranslationBounds OffsetBoundsToTranslationBounds.cc)
target_link_libraries(OffsetBoundsToTranslationBounds isosurface_inflator pattern_optimization)

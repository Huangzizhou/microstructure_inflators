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
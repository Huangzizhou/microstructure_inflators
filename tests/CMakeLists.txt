################################################################################
cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(unit_tests
	main.cc
	test_isosurfacejob.cc
	test_inflation.cc
	test_stitching.cc
)

target_link_libraries(unit_tests PUBLIC
	Catch2::Catch2
	isosurface_inflator
	pattern_optimization
)

# Data directory
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data/")
target_compile_definitions(unit_tests PUBLIC -DDATA_DIR=\"${DATA_DIR}\")

# Register tests
include(CTest)
include(Catch)
catch_discover_tests(unit_tests)

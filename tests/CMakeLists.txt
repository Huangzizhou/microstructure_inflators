################################################################################
cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(unit_tests
	main.cpp
	test_stitching.cpp
)

target_link_libraries(unit_tests PUBLIC
	Catch2::Catch2
	isosurface_inflator
)

# Register tests
include(CTest)
include(Catch)
catch_discover_tests(unit_tests)

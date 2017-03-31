################################################################################
# General Informations
################################################################################

cmake_minimum_required(VERSION 3.1)
project(Microstructures)

# CMP0063: Honor visibility properties for all target types.
if(POLICY CMP0063)
	cmake_policy(SET CMP0063 NEW)
endif()

################################################################################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Color output
include(UseColors)

# Export compile flags(used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# CMake plugin for vscode
include(CMakeToolsHelpers OPTIONAL)

################################################################################
# Sub-projects
################################################################################

add_subdirectory(isosurface_inflator)
add_subdirectory(pattern_optimization)
add_subdirectory(worst_case_stress)

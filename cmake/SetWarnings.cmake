################################################################################
cmake_minimum_required(VERSION 2.6.3)
################################################################################
# See comments and discussions here:
# http://stackoverflow.com/questions/5088460/flags-to-enable-thorough-and-verbose-g-warnings
################################################################################

set(MY_FLAGS
		-Wall
		-Wextra
		-Wpedantic
		-Wno-comment

		# Gives meaningful stack traces
		-fno-omit-frame-pointer
)

# Flags above don't make sense for MSVC
if(MSVC)
	set(MY_FLAGS)
endif()

include(CheckCXXCompilerFlag)

set(ALL_WARNINGS)
foreach(FLAG ${MY_FLAGS})
	check_cxx_compiler_flag("${FLAG}" IS_SUPPORTED_${FLAG})
	if(IS_SUPPORTED_${FLAG})
		set(ALL_WARNINGS ${ALL_WARNINGS} ${FLAG})
	endif()
endforeach()

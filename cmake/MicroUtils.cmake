# Create single-source application
function(micro_single_app name)
	add_executable(${name} ${name}.cc)
	target_link_libraries(${name} ${ARGN})
endfunction()

# Copy header files into the target folder of the build directory
# Only works with relatives paths so far
function(micro_copy_headers target)
	if(NOT TARGET ${target})
		message(WARNING "${target} is not a CMake target. micro_copy_headers() will not be run.")
		return()
	endif()

	# Get info on target
	get_target_property(TARGET_NAME ${target} NAME)
	get_target_property(TARGET_TYPE ${target} TYPE)
	if(TARGET_TYPE STREQUAL INTERFACE_LIBRARY)
		set(TARGET_SCOPE INTERFACE)
	else()
		set(TARGET_SCOPE PUBLIC)
	endif()

	if(MICRO_COPY_HEADERS)
		# Copy header files
		foreach(filepath IN ITEMS ${ARGN})
			set(filename "${filepath}")
			# get_filename_component(filename "${filepath}" NAME)
			if(${filename} MATCHES ".*\.(hh|h|inl)$")
				configure_file(${filepath} ${PROJECT_BINARY_DIR}/include/${TARGET_NAME}/${TARGET_NAME}/${filename})
			endif()
		endforeach()

		# Set target include directory
		target_include_directories(${target} ${TARGET_SCOPE} ${PROJECT_BINARY_DIR}/include/${TARGET_NAME})
	else()
		# Set target include directory
		target_include_directories(${target} ${TARGET_SCOPE} ..)
	endif()

endfunction()

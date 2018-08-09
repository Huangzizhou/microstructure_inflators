################################################################################
# Find PyMesh sources
# The following are set:
#
# PYMESH_FOUND     - Whether the PyMesh and PyWires libraries were found
# PyMesh::core     - Imported target for PyMesh
# PyMesh::wires    - Imported target for the Wires library
#
# It searches the environment variable $PYMESH_PATH
################################################################################

set(_PYMESH_SEARCH_PATHS
		"$ENV{PYMESH_PATH}"
		${MICRO_EXTERNAL}/PyMesh
		"C:/Program Files/PyMesh/"
		"$ENV{HOME}/external/git/PyMesh"
)

find_path(PYMESH_INCLUDE
		"Geometry/MeshGeometry.h"
		PATHS
			${_PYMESH_SEARCH_PATHS}
		PATH_SUFFIXES
			src
)

find_library(PYMESH_LIBRARY
		NAMES Mesh
		PATHS ${_PYMESH_SEARCH_PATHS}
		PATH_SUFFIXES python/pymesh/lib
)

find_library(PYWIRES_LIBRARY
		NAMES wires
		PATHS ${_PYMESH_SEARCH_PATHS}
		PATH_SUFFIXES python/pymesh/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PYMESH DEFAULT_MSG PYMESH_INCLUDE PYMESH_LIBRARY PYWIRES_LIBRARY)

if(PYMESH_FOUND AND NOT TARGET PyMesh::core)
	# Imported interface target for PyMesh
	add_library(PyMesh::core UNKNOWN IMPORTED)

	# Interface include directory
	set_target_properties(PyMesh::core PROPERTIES
			INTERFACE_INCLUDE_DIRECTORIES "${PYMESH_INCLUDE}")

	# Link to library file
	set_target_properties(PyMesh::core PROPERTIES
			IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
			IMPORTED_LOCATION "${PYMESH_LIBRARY}")
endif()

if(PYMESH_FOUND AND NOT TARGET PyMesh::wires)
	# Imported interface target for PyMesh
	add_library(PyMesh::wires UNKNOWN IMPORTED)

	# Interface include directory
	set_target_properties(PyMesh::wires PROPERTIES
			INTERFACE_INCLUDE_DIRECTORIES "${PYMESH_INCLUDE}/../tools")

	# Link to library file
	set_target_properties(PyMesh::wires PROPERTIES
			IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
			IMPORTED_LOCATION "${PYWIRES_LIBRARY}")
endif()

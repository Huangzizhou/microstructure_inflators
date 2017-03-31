# - Find Clipper
# Find the clipper includes and library
#  CLIPPER_INCLUDE_DIR - where to find clipper.hpp
#  CLIPPER_LIBRARIES   - clipper library binary
#  CLIPPER_FOUND       - True if clipper found.


IF (CLIPPER_INCLUDE_DIR)
  # Already in cache, be silent
  SET (CLIPPER_FIND_QUIETLY TRUE)
ENDIF (CLIPPER_INCLUDE_DIR)

FIND_PATH(CLIPPER_INCLUDE_DIR clipper.hpp
	PATHS
		"C:/Program Files/Clipper/include"
		$ENV{CLIPPER_INC}
)

FIND_LIBRARY (CLIPPER_LIBRARIES polyclipping
	PATHS
		"C:/Program Files/Clipper/lib"
		$ENV{CLIPPER_LIB}
)

# handle the QUIETLY and REQUIRED arguments and set CLIPPER_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (CLIPPER DEFAULT_MSG
  CLIPPER_LIBRARIES
  CLIPPER_INCLUDE_DIR)

MARK_AS_ADVANCED (CLIPPER_LIBRARIES CLIPPER_INCLUDE_DIR)

#
# The following module is based on FindVTK.cmake
#

# - Find a CGAL installation or binary tree.
# The following variables are set if CGAL is found.  If CGAL is not
# found, CGAL_FOUND is set to false.
#
#  CGAL_FOUND         - Set to true when CGAL is found.
#

# Construct consistent error messages for use below.
set(CGAL_DIR_DESCRIPTION "directory containing CGALConfig.cmake. This is either the binary directory where CGAL was configured or PREFIX/lib/CGAL for an installation.")
set(CGAL_DIR_MESSAGE     "CGAL not found. Set the CGAL_DIR cmake variable or environment variable to the ${CGAL_DIR_DESCRIPTION}")

include(CMakeFindDependencyMacro)
find_dependency(CGAL
    CONFIG
    PATHS
        # Look for an environment variable CGAL_DIR.
        $ENV{CGAL_DIR}
        $ENV{CGAL_ROOT}

        # Look in the standard Macports install location
        /opt/local/share/CGAL/cmake

    # Help the user find it if we cannot.
    PATH_SUFFIXES lib64/CGAL
)

################################################################################
# Find MeshFEM
# The following are set:
#
# MeshFEM_FOUND - Whether the MeshFEM library was found
# MeshFEM       - Target to build MeshFEM along with the project
#
# It searches the environment variable $MESHFEM_PATH
################################################################################

find_path(MESHFEM_INCLUDE
		FEMMesh.hh
		HINTS
			${PROJECT_SOURCE_DIR}/../MeshFEM
			${PROJECT_SOURCE_DIR}/../../jpanetta/MeshFEM
		PATHS
			ENV MESHFEM_PATH
			${THIRD_PARTY_DIR}/MeshFEM
			"C:/Program Files/MeshFEM/"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MESHFEM DEFAULT_MSG MESHFEM_INCLUDE)

if(MESHFEM_FOUND AND NOT TARGET MeshFEM)
	# Build MeshFEM library alongside the rest of the project
	add_subdirectory(${MESHFEM_INCLUDE} MeshFEM)
endif()

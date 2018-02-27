# - Find a Ceres installation.

if(NOT Ceres_DIR)

  find_path(Ceres_DIR CeresConfig.cmake
    PATHS
      ${THIRD_PARTY_DIR}/ceres
    PATH_SUFFIXES
      lib/cmake/Ceres
      lib64/cmake/Ceres
  )

endif()

find_package(Ceres CONFIG)

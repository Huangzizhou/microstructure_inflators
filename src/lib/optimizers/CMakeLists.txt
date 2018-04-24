################################################################################

add_library(optimizers
        wrappers/ceres.cc
        wrappers/dlib.cc
        wrappers/gradient_descent.cc
        wrappers/nlopt.cc
        wrappers/knitro.cc
)

################################################################################
# Required libraries
################################################################################

# MeshFEM + inflators
target_link_libraries(optimizers MeshFEM inflators)

# Accelerate framework
if(APPLE)
    find_library(AccelerateFramework Accelerate)
    target_link_libraries(optimizers ${AccelerateFramework})
endif()

################################################################################
# Optional libraries
################################################################################

# Dlib Library
find_package(Dlib QUIET)
if(DLIB_FOUND)
    target_include_directories(optimizers SYSTEM PUBLIC ${DLIB_INCLUDE_DIR})
    target_link_libraries(optimizers ${DLIB_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_DLIB)
else()
    message(STATUS "DLib not found; bfgs optimizer is disabled")
endif()

# Ceres library
find_package(Ceres QUIET)
if(CERES_FOUND)
    target_include_directories(optimizers SYSTEM PUBLIC ${CERES_INCLUDE_DIRS})
    target_link_libraries(optimizers ${CERES_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_CERES)
else()
    message(STATUS "Google's ceres-solver not found; levenberg-marquardt disabled")
endif()

# NLopt library
find_package(NLopt QUIET)
if(NLOPT_FOUND)
    target_link_libraries(optimizers ${NLOPT_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_NLOPT)
else()
    message(STATUS "NLopt not found; slsqp disabled")
endif()

# Knitro library
find_package(Knitro QUIET)
if(KNITRO_FOUND)
    target_include_directories(optimizers SYSTEM PUBLIC ${KNITRO_INCLUDE_DIRS})
    target_link_libraries(optimizers ${KNITRO_LIBRARIES})
    target_compile_definitions(optimizers PUBLIC -DHAS_KNITRO)
else()
    message(STATUS "Knitro not found; active_set disabled")
endif()

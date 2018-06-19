################################################################################
# Source files
################################################################################

add_library(optimizers
        wrappers/ceres.cc
        wrappers/dlib.cc
        wrappers/gradient_descent.cc
        wrappers/nlopt.cc
        wrappers/knitro.cc
)

################################################################################
# Dependencies
################################################################################

target_link_libraries(optimizers
        MeshFEM
        inflators
        micro::accelerate
        micro::dlib
        micro::ceres
        micro::nlopt
        micro::knitro
)

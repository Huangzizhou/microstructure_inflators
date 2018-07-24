################################################################################
# Source files
################################################################################

set(SOURCES
        BoundConstraints.hh
        Constraint.hh
        EvaluatedConstraint.hh
        EvaluatedObjectiveTerm.hh
        IterateBase.hh
        IterateManagerBase.hh
        OptimizerConfig.hh
        wrappers/ceres.cc
        wrappers/ceres.hh
        wrappers/dlib.cc
        wrappers/dlib.hh
        wrappers/gradient_descent.cc
        wrappers/gradient_descent.hh
        wrappers/knitro.cc
        wrappers/knitro.hh
        wrappers/nlopt.cc
        wrappers/nlopt.hh
)

add_library(optimizers ${SOURCES})
micro_copy_headers(optimizers ${SOURCES})

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

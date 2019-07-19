################################################################################
# Source files
################################################################################

set(SOURCES
        MicroscopicStress.hh
        NonPeriodicCellOperations.hh
        ParametersMask.hh
        PolygonMesher.cc
        PolygonMesher.hh
        ShapeOptimizationIterate.hh
        StressObjectiveTerm.hh
)

micro_add_library(shape_optimization ${SOURCES})
micro_copy_headers(shape_optimization ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_shape_optimization
    PUBLIC
        MeshFEM
        meshfem::boost
        micro::isosurface_inflator
    PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(micro_shape_optimization)
endif()

################################################################################
# Source files
################################################################################

set(SOURCES
        Algebra.hh
        BaseCellType.hh
        BoundaryMesh.hh
        CollisionGrid.hh
        ComponentMask.hh
        Concepts.hh
        Fields.hh
        Flattening.hh
        function_traits.hh
        Functions.hh
        Future.hh
        Geometry.hh
        GlobalBenchmark.cc
        GlobalBenchmark.hh
        MeshDataTraits.hh
        MeshIO.cc
        MeshIO.hh
        MSHFieldParser.cc
        MSHFieldParser.hh
        MSHFieldWriter.hh
        NTuple.hh
        Parallelism.hh
        PeriodicBoundaryMatcher.hh
        Simplex.hh
        SimplicialMesh.hh
        SimplicialMeshInterface.hh
        StringUtils.cc
        StringUtils.hh
        SymmetricMatrix.hh
        TemplateHacks.hh
        TetMesh.hh
        TetMesh.inl
        Timer.hh
        Triangulate.h
        TriMesh.hh
        TriMesh.inl
        Types.cc
        Types.hh
        util.h
        utils.hh

        algorithms/get_element_components.hh
        algorithms/remove_if_index.hh
        filters/CurveCleanup.hh
        filters/extract_polygons.hh
        filters/gen_grid.hh
        filters/reflect.hh
        filters/remove_dangling_vertices.hh
        filters/remove_small_components.hh
        filters/ResampleCurve.hh
        Handles/Handle.hh
        Handles/TetMeshHandles.hh
        Handles/TriMeshHandles.hh
        Utilities/apply.hh
        Utilities/EdgeAccessAdaptor.hh
        Utilities/EdgeSoupAdaptor.hh
        Utilities/IteratorMap.hh
        Utilities/NDArray.hh
        Utilities/RandomAccessIndexSet.hh
        wrappers/meshfem_triangle.h
)

micro_add_library(mesh_fem ${SOURCES})
micro_copy_headers(mesh_fem ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_mesh_fem
    PUBLIC
        Threads::Threads
        Eigen3::Eigen
        json::json
        micro::tbb
        optional::optional
        triangle::triangle
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(MeshFEM)
endif()

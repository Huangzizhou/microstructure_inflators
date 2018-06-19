#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import glob

file_list = [
    "MeshFEM/Geometry.hh",
    "MeshFEM/Simplex.hh",
    "MeshFEM/BaseCellType.hh",
    "MeshFEM/MaterialField.hh",
    "MeshFEM/OneForm.hh",
    "MeshFEM/MeshIO.hh",
    "MeshFEM/MassMatrix.hh",
    "MeshFEM/MSHFieldWriter.hh",
    "MeshFEM/NTuple.hh",
    "MeshFEM/Laplacian.hh",
    "MeshFEM/InterpolantRestriction.hh",
    "MeshFEM/ExpressionVector.hh",
    "MeshFEM/GlobalBenchmark.hh",
    "MeshFEM/Concepts.hh",
    "MeshFEM/BoundaryMesh.hh",
    "MeshFEM/SymmetricMatrixInterpolant.hh",
    "MeshFEM/GridFunction.hh",
    "MeshFEM/TetMesh.hh",
    "MeshFEM/DenseCollisionGrid.hh",
    "MeshFEM/MaterialOptimization.hh",
    "MeshFEM/CollisionGrid.hh",
    "MeshFEM/LinearElasticity.hh",
    "MeshFEM/Triangulate.h",
    "MeshFEM/MeshDataTraits.hh",
    "MeshFEM/FEMMesh.hh",
    "MeshFEM/PerturbMesh.hh",
    "MeshFEM/Functions.hh",
    "MeshFEM/PeriodicHomogenization.hh",
    "MeshFEM/MSHFieldParser.hh",
    "MeshFEM/BoundaryConditions.hh",
    "MeshFEM/Future.hh",
    "MeshFEM/OrthotropicHomogenization.hh",
    "MeshFEM/TriMesh.hh",
    "MeshFEM/Types.hh",
    "MeshFEM/HalfedgeDictionary.hh",
    "MeshFEM/PeriodicBoundaryMatcher.hh",
    "MeshFEM/GaussQuadrature.hh",
    "MeshFEM/TemplateHacks.hh",
    "MeshFEM/filters/extrude.hh",
    "MeshFEM/filters/quad_tri_subdiv_asymmetric.hh",
    "MeshFEM/filters/quad_tri_subdiv.hh",
    "MeshFEM/filters/ResampleCurve.hh",
    "MeshFEM/filters/quad_subdiv_high_aspect.hh",
    "MeshFEM/filters/hex_tet_subdiv.hh",
    "MeshFEM/filters/extract_polygons.hh",
    "MeshFEM/filters/gen_cursor.hh",
    "MeshFEM/filters/CurveCleanup.hh",
    "MeshFEM/filters/subdivide.hh",
    "MeshFEM/filters/extract_hole_boundaries.hh",
    "MeshFEM/filters/remove_dangling_vertices.hh",
    "MeshFEM/filters/gen_grid.hh",
    "MeshFEM/filters/remove_small_components.hh",
    "MeshFEM/filters/quad_subdiv.hh",
    "MeshFEM/filters/reorient_negative_elements.hh",
    "MeshFEM/filters/merge_duplicate_vertices.hh",
    "MeshFEM/filters/reflect.hh",
    "MeshFEM/filters/voxels_to_simplices.hh",
    "MeshFEM/filters/highlight_dangling_vertices.hh",
    "MeshFEM/JSFieldWriter.hh",
    "MeshFEM/algorithms/remove_if_index.hh",
    "MeshFEM/algorithms/get_element_components.hh",
    "MeshFEM/Poisson.hh",
    "MeshFEM/EmbeddedElement.hh",
    "MeshFEM/HalfEdge.hh",
    "MeshFEM/Utilities/NDArray.hh",
    "MeshFEM/Utilities/IteratorMap.hh",
    "MeshFEM/Utilities/EdgeSoupAdaptor.hh",
    "MeshFEM/Utilities/EdgeAccessAdaptor.hh",
    "MeshFEM/Utilities/ci_string.hh",
    "MeshFEM/Utilities/RandomAccessIndexSet.hh",
    "MeshFEM/Utilities/apply.hh",
    "MeshFEM/EdgeFields.hh",
    "MeshFEM/SimplicialMesh.hh",
    "MeshFEM/ComponentMask.hh",
    "MeshFEM/BoundaryLaplacian.hh",
    "MeshFEM/Materials.hh",
    "MeshFEM/function_traits.hh",
    "MeshFEM/SimplicialMeshInterface.hh",
    "MeshFEM/util.h",
    "MeshFEM/UniformLaplacian.hh",
    "MeshFEM/Handles/FEMMeshHandles.hh",
    "MeshFEM/Handles/Handle.hh",
    "MeshFEM/Handles/TriMeshHandles.hh",
    "MeshFEM/Handles/TetMeshHandles.hh",
    "CSGFEM/CSGWindowController.hh",
    "CSGFEM/ShaderCompiler.hh",
    "CSGFEM/TensorProjection.hh",
    "CSGFEM/Geometry.hh",
    "CSGFEM/Flipbook.hh",
    "CSGFEM/ParameterSweepDialog.hh",
    "CSGFEM/Solver.hh",
    "CSGFEM/CSGWindow.hh",
    "CSGFEM/SparseMatrices.hh",
    "CSGFEM/SymmetricMatrix.hh",
    "CSGFEM/Fields.hh",
    "CSGFEM/CSGTreeModel.hh",
    "CSGFEM/Timer.hh",
    "CSGFEM/LinearIndexer.hh",
    "CSGFEM/ElementGrid.hh",
    "CSGFEM/Algebra.hh",
    "CSGFEM/Flattening.hh",
    "CSGFEM/precompile.hh",
    "CSGFEM/colors.hh",
    "CSGFEM/SPHKernels.hh",
    "CSGFEM/ViewSettings.hh",
    "CSGFEM/CSGFile.hh",
    "CSGFEM/ViewSettingsWidget.hh",
    "CSGFEM/BoundaryConditionDialog.hh",
    "CSGFEM/LazyMatlabInterfaces.hh",
    "CSGFEM/Quadrature.hh",
    "CSGFEM/ParameterSweep.hh",
    "CSGFEM/BoundaryConditions.hh",
    "CSGFEM/build_specs/macx-g++47/qplatformdefs.h",
    "CSGFEM/GlobalTypes.hh",
    "CSGFEM/SolverLibrary.hh",
    "CSGFEM/AnalysisSettings.hh",
    "CSGFEM/VonMises.hh",
    "CSGFEM/ModelForm.hh",
    "CSGFEM/LevelSet.hh",
    "CSGFEM/QuadraturePointsSpinBox.hh",
    "CSGFEM/Parallelism.hh",
    "CSGFEM/draw.hh",
    "CSGFEM/MeshlessFEM2D.hh",
    "CSGFEM/AnalysisForm.hh",
    "CSGFEM/QCommandLine.hh",
    "CSGFEM/ResultsCollector.hh",
    "CSGFEM/QMatlabInterface.hh",
    "CSGFEM/MatlabInterface/MatlabShell.h",
    "CSGFEM/MatlabInterface/MatlabInterface.h",
    "CSGFEM/WireNetwork.hh",
    "CSGFEM/FEMView.hh",
    "CSGFEM/Grid.hh",
    "CSGFEM/MarchingSquaresGrid.hh",
    "CSGFEM/CSGTree.hh",
    "CSGFEM/MSHWriter.hh",
    "CSGFEM/BenchmarkStub.hh",
    "CSGFEM/MeshlessFEM3D.hh",
    "CSGFEM/ElasticityTensor.hh",
    "CSGFEM/cl-helper.h",
    "CSGFEM/ResultsWindow/ResultTreeView.hh",
    "CSGFEM/ResultsWindow/ResultsWindow.hh",
    "CSGFEM/ResultsWindow/ResultsWindowController.hh",
    "CSGFEM/ResultsWindow/FlipbookDialog.hh",
    "CSGFEM/utils.hh",
]

all_keys = { x.partition('/')[2] for x in file_list }


def find_match(key):
    for x in file_list:
        if x.endswith(key):
            return x
    return key


def process_file(filename):
    print(filename)
    with open(filename, 'r') as f:
        content = f.read()
    tokens = re.findall('^#include <(.*)>', content, re.MULTILINE)
    for key in tokens:
        if key in all_keys:
            str_from = '#include <{}>'.format(key)
            str_to = '#include <{}>'.format(find_match(key))
            content = content.replace(str_from, str_to)
    with open(filename, 'w') as f:
        f.write(content)


def process_folder(folder):
    types = ('*.c', '*.cc', '*.cpp', '*.h', '*.hh', '*.hpp', '*.inl')
    files_grabbed = []
    for pattern in types:
        files_grabbed.extend(glob.glob(os.path.join(folder, pattern)))
    for file in files_grabbed:
        process_file(file)


def process_folders(folder):
    for subdir, _dirs, files in os.walk(folder):
        process_folder(subdir)

if __name__ == "__main__":
    pass
    # process_folders(folder)

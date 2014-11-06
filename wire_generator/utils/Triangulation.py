import numpy as np
from numpy.linalg import norm

import PyMesh
import PyMeshUtils
import PyTriangle

class Triangulation(object):
    def __init__(self, vertices, faces, project_dir):
        assert(vertices.ndim == 2 and vertices.shape[1] == 3);
        assert(faces.ndim == 2 and faces.shape[1] == 3);

        self.ori_vertices = vertices;
        self.ori_faces = faces;
        self.project_dir = project_dir / norm(project_dir);

        self.__extract_boundary();
        self.__remove_isolated_bd_nodes();

    def triangulate(self, max_area=0.1):
        tri = PyTriangle.TriangleWrapper(self.bd_nodes, self.bd_edges);
        tri.run(max_area, False, True);

        self.vertices = tri.get_vertices();
        self.faces = tri.get_faces();

    def __extract_boundary(self):
        bd_extractor = PyMeshUtils.Boundary.extract_surface_boundary_raw(
                self.ori_vertices, self.ori_faces);
        self.bd_nodes = self.ori_vertices;
        self.bd_edges = bd_extractor.get_boundaries();

    def __remove_isolated_bd_nodes(self):
        remover = PyMeshUtils.IsolatedVertexRemoval(
                self.bd_nodes, self.bd_edges);
        remover.run();

        self.bd_nodes = remover.get_vertices();
        self.bd_edges = remover.get_faces();

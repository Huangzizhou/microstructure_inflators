import numpy as np
from numpy.linalg import norm

import PyMesh
import PyMeshUtils
import triangle

class Triangulation(object):
    def __init__(self, vertices, faces, project_dir):
        assert(vertices.ndim == 2 and vertices.shape[1] == 3);
        assert(faces.ndim == 2 and faces.shape[1] == 3);

        self.ori_vertices = vertices;
        self.ori_faces = faces;
        self.project_dir = project_dir / norm(project_dir);

        self.__compute_2D_coordinate_frame();
        self.__remove_isolated_vertices();
        self.__extract_boundary_edge_loops();

    def __remove_isolated_vertices(self):
        remover = PyMeshUtils.IsolatedVertexRemoval(
                self.ori_vertices, self.ori_faces);
        remover.run();
        self.ori_vertices = remover.get_vertices();
        self.ori_faces = remover.get_faces();

    def __extract_boundary_edge_loops(self):
        mesh = self._form_mesh(self.ori_vertices, self.ori_faces);
        bd_extractor = PyMeshUtils.Boundary.extract_surface_boundary(mesh);

        self.bd_nodes = self.ori_vertices[bd_extractor.get_boundary_nodes().ravel()];
        self.bd_edges = bd_extractor.get_boundaries();

        #self.bd_loops = self.__extract_edge_loops(self.bd_nodes, self.bd_edges);

    def __project_onto_2D_plane(self, vertices):
        bd_nodes_coords = np.dot(self.coordinate_frame, vertices.T).T;
        return bd_nodes_coords[:,:2], np.average(bd_nodes_coords[:,2]);

    def __project_onto_3D_plane(self, vertices, intercept):
        intercept = np.ones((len(vertices), 1)) * intercept;
        vertices = np.hstack([vertices, intercept]);
        return np.dot(self.coordinate_frame.T, vertices.T).T;

    def triangulate(self, max_area=0.1):
        self.proj_bd_nodes, intercept =\
                self.__project_onto_2D_plane(self.bd_nodes);
        tri_mesh = {
                "vertices": self.proj_bd_nodes,
                "segments": self.bd_edges,
                }
        tri_options = "pqQYa{}".format(max_area);
        result = triangle.triangulate(tri_mesh, tri_options);

        self.vertices = result["vertices"];
        self.vertices = self.__project_onto_3D_plane(self.vertices, intercept);
        self.faces = result["triangles"];

    def __compute_2D_coordinate_frame(self):
        cononical_frame = np.eye(3);

        for axis_dir in cononical_frame:
            proj_len = np.dot(axis_dir, self.project_dir);
            if abs(proj_len) > 1.0 - 1e-3: continue;

            z_dir = self.project_dir;
            x_dir = axis_dir - self.project_dir * proj_len;
            x_dir /= norm(x_dir);
            y_dir = np.cross(z_dir, x_dir);
            self.coordinate_frame = np.array([x_dir, y_dir, z_dir]);
            return;


    def __extract_edge_loops(self, vertices, edges):
        num_vertices = len(vertices);
        adj_list = np.ones((num_vertices, 2), dtype=int) * -1;
        for edge in edges:
            if adj_list[edge[0], 0] < 0:
                adj_list[edge[0], 0] = edge[1];
            elif adj_list[edge[0], 1] < 0:
                adj_list[edge[0], 1] = edge[1];
            else:
                raise NotImplementedError("Non-manifold loop detected");

            if adj_list[edge[1], 0] < 0:
                adj_list[edge[1], 0] = edge[0];
            elif adj_list[edge[1], 1] < 0:
                adj_list[edge[1], 1] = edge[0];
            else:
                raise NotImplementedError("Non-manifold loop detected");

        loops = [];
        visited = np.zeros(num_vertices, dtype=bool);
        for i in range(num_vertices):
            if visited[i]: continue;
            loop = [i];
            assert(not visited[adj_list[i, 0]]);
            loop.append(adj_list[i, 0]);

            while True:
                neighbor = adj_list[loop[-1]];
                if neighbor[0] == loop[-2]:
                    assert(not visited[neighbor[1]]);
                    next_v = neighbor[1];
                else:
                    assert(not visited[neighbor[0]]);
                    next_v = neighbor[0];
                if loop[0] == next_v:
                    break;
                else:
                    loop.append(next_v);

            visited[loop] = True;
            loops.append(loop);
        return loops;

    def _form_mesh(self, vertices, faces):
        voxels = np.array([]);
        factory = PyMesh.MeshFactory();
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels,
                3, 3, 0);
        return factory.create();


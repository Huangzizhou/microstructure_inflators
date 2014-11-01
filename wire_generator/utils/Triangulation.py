import numpy as np
from numpy.linalg import norm, inv
from math import pi

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

        #self.bd_nodes = self.ori_vertices[bd_extractor.get_boundary_nodes().ravel()];
        self.bd_nodes = self.ori_vertices;
        self.bd_edges = bd_extractor.get_boundaries();

        self.__extract_edge_loops();

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
        hole_points = self.__compute_hole_points().reshape((-1, 2));
        tri_mesh = {
                "vertices": self.proj_bd_nodes,
                "segments": self.bd_edges,
                }
        if len(hole_points) > 0:
            tri_mesh["holes"] = hole_points;
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

    def __compute_hole_points(self):
        self.__find_interior_points();
        self.__compute_loop_orientations();
        return self.interior_points[np.logical_not(self.loop_orientations)];

    def __compute_adjacency_list(self):
        vertices = self.bd_nodes;
        edges = self.bd_edges;
        num_vertices = len(vertices);
        adj_list = np.ones((num_vertices, 2), dtype=int) * -1;
        for edge in edges:
            assert(adj_list[edge[0], 0] < 0);
            assert(adj_list[edge[1], 1] < 0);
            adj_list[edge[0], 0] = edge[1];
            adj_list[edge[1], 1] = edge[0];
        return adj_list;

    def __extract_edge_loops(self):
        vertices = self.bd_nodes;
        edges = self.bd_edges;
        num_vertices = len(vertices);

        adj_list = self.__compute_adjacency_list();

        loops = [];
        visited = np.zeros(num_vertices, dtype=bool);
        for i in range(num_vertices):
            if visited[i]: continue;
            loop = [i];
            assert(not visited[adj_list[i, 0]]);
            loop.append(adj_list[i, 0]);

            while True:
                neighbor = adj_list[loop[-1]];
                next_v = neighbor[0];
                assert(not visited[next_v]);
                assert(loop[-2] != next_v);

                if loop[0] == next_v:
                    break;
                else:
                    loop.append(next_v);

            visited[loop] = True;
            loops.append(loop);
        self.bd_loops = loops;

    def __find_interior_points(self):
        self.interior_points = np.array([self.__find_interior_point(loop)
            for loop in self.bd_loops]);

    def __compute_loop_orientations(self):
        self.loop_orientations = [];
        for i,loop in enumerate(self.bd_loops):
            vertices = self.proj_bd_nodes[loop]
            p = self.interior_points[i];
            wind_num = self.__compute_winding_number(p, vertices);
            self.loop_orientations.append(wind_num > 0.0);

    def __compute_winding_number(self, p, vertices):
        num_vts = len(vertices);
        z_coord = np.zeros((num_vts,1));
        v0 = np.hstack((vertices - p, z_coord));
        v1 = np.roll(v0, -1, axis=0);
        angles = np.arctan2(
                np.cross(v0, v1)[:,2],
                np.sum(v0 * v1, axis=1));
        assert(len(angles) == num_vts);
        wind_num = np.sum(angles) / (2*pi);
        return wind_num;

    def __find_interior_point(self, loop):
        """ For a closed loop, find a point that is interior. Based on
        modification of the following algorithm:
        http://www.alecjacobson.com/weblog/?p=1256
        """
        vertices = self.proj_bd_nodes[loop];
        bbox_min = np.amin(vertices, axis=0);
        bbox_max = np.amax(vertices, axis=0);
        center = 0.5 * (bbox_min + bbox_max);

        rotation = np.array([[0, 1], [-1, 0]]);

        t0 = center + (bbox_max - center) * 2;
        t1 = center - (bbox_max - center) * 2;
        t_normal = np.dot(rotation, t1 - t0);

        p0 = vertices;
        p1 = np.roll(vertices, -1, axis=0);
        p_normal = np.dot(rotation, (p1 - p0).T).T;

        t_sides = [
                np.dot(t_normal, (p0 - t0).T).T,
                np.dot(t_normal, (p1 - t0).T).T ];
        p_sides = [
                np.sum(p_normal * (t0 - p0), axis=1),
                np.sum(p_normal * (t1 - p0), axis=1) ];
        crossed = np.logical_and(
                np.logical_xor(t_sides[0] < 0, t_sides[1] < 0),
                np.logical_xor(p_sides[0] < 0, p_sides[1] < 0));
        assert(np.sum(crossed) >= 2);

        p0 = p0[crossed];
        p1 = p1[crossed];
        crossed_points = [];
        for q0, q1 in zip(p0, p1):
            sol = np.dot(inv(np.array([t1-t0, q0-q1]).T), (q0 - t0));
            p = t0 + sol[0] * (t1 - t0);
            crossed_points.append(p);
        crossed_points = np.array(crossed_points);

        dist_to_t0 = norm(crossed_points - t0, axis=1);
        order = np.argsort(dist_to_t0);
        return 0.5 * (crossed_points[order[0]] + crossed_points[order[1]]);

    def _form_mesh(self, vertices, faces):
        voxels = np.array([]);
        factory = PyMesh.MeshFactory();
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels,
                3, 3, 4);
        return factory.create();


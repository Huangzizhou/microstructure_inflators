import os
import sys

import numpy as np
from numpy.linalg import norm
from math import pi, sqrt
from scipy.spatial import ConvexHull

from WireNetwork import WireNetwork
from WireInflator import WireInflator
from WirePattern import WirePattern
from timethis import timethis

import PyMesh
import PyMeshUtils

class PeriodicWireInflator(WireInflator):
    def __init__(self, wire_network):
        super(PeriodicWireInflator, self).__init__(wire_network);
        self.__generate_phantom_wire_network();
        self.__generate_phantom_vertex_map();

        # Use phantom wires as input for inflator
        self.original_wire_network = self.wire_network;
        self.wire_network = self.phantom_wire_network;


    def inflate(self, thickness, clean_up=True):
        if not clean_up:
            raise NotImplementedError("Clean up is required for periodic wires");
        super(PeriodicWireInflator, self).inflate(thickness, clean_up);
        self._clip_exterior();
        self._enforce_periodic_connectivity();

    def __generate_phantom_wire_network(self):
        wire_pattern = WirePattern();
        wire_pattern.set_single_cell_from_wire_network(self.wire_network);
        wire_pattern.tile([3, 3, 3]);
        phantom_wire_network = wire_pattern.wire_network;

        phantom_wire_network.translate(
                self.wire_network.centroid - phantom_wire_network.centroid);
        self.phantom_wire_network = phantom_wire_network;

    def __generate_phantom_vertex_map(self):
        bbox_min, bbox_max = self.phantom_wire_network.bbox;
        cell_size = np.amin(bbox_max - bbox_min) * 0.001;
        hash_grid = PyMesh.HashGrid.create(cell_size);
        hash_grid.insert_multiple(
                np.arange(len(self.wire_network.vertices), dtype=int),
                self.wire_network.vertices);

        self.phantom_vertex_map = -1 * np.ones(
                self.phantom_wire_network.num_vertices, dtype=int);
        for i,v in enumerate(self.phantom_wire_network.vertices):
            nearby_v_indices = hash_grid.get_items_near_point(v);
            if len(nearby_v_indices) > 0:
                if len(nearby_v_indices) > 1:
                    raise RuntimeError(
                            "Periodic unit pattern contains duplicated vertices");
                self.phantom_vertex_map[i] = nearby_v_indices[0];

    def _generate_edge_pipes(self):
        for i,edge in enumerate(self.wire_network.edges):
            if np.any(self.phantom_vertex_map[edge] < 0): continue;
            loop_1_idx = self.edge_loop_indices[2*i  , :];
            loop_2_idx = self.edge_loop_indices[2*i+1, :];
            faces = [
                    [loop_1_idx[0], loop_1_idx[1], loop_2_idx[0]],
                    [loop_2_idx[0], loop_1_idx[1], loop_2_idx[1]],
                    [loop_2_idx[1], loop_1_idx[1], loop_1_idx[2]],
                    [loop_2_idx[1], loop_1_idx[2], loop_2_idx[2]],
                    [loop_2_idx[2], loop_1_idx[2], loop_1_idx[3]],
                    [loop_2_idx[2], loop_1_idx[3], loop_2_idx[3]],
                    [loop_2_idx[3], loop_1_idx[3], loop_1_idx[0]],
                    [loop_2_idx[3], loop_1_idx[0], loop_2_idx[0]] ];
            self.mesh_faces = np.vstack((self.mesh_faces, faces));

    def _generate_joint(self, idx):
        if self.phantom_vertex_map[idx] < 0:
            return;

        num_vts = len(self.mesh_vertices);
        v = self.wire_network.vertices[idx];
        loop_vertices = v.reshape((-1,3));
        phantom_loop = [];
        for i,e_idx in enumerate(self.wire_network.vertex_edge_neighbors[idx]):
            edge = self.wire_network.edges[e_idx];
            v_idx = np.where(np.array(edge).ravel() == idx)[0][0];
            loop = self.edge_loops[e_idx, v_idx, :, :].reshape((-1,3), order="C");
            loop_vertices = np.vstack((loop_vertices, loop));
            self.edge_loop_indices[e_idx * 2 + v_idx, :] = \
                    np.arange(i*4, i*4+4) + num_vts + 1;
            phantom_loop.append(np.any(self.phantom_vertex_map[edge] < 0));

        def is_edge_loop_face(simplex):
            idx = [(vi - 1) / 4 for vi in simplex];
            if np.max(idx) == np.min(idx):
                # all face vertices belong to the same edge loop
                if not phantom_loop[idx[0]]:
                    return True;
            return False;

        joint_center = np.mean(loop_vertices, axis=0);
        convex_hull = ConvexHull(loop_vertices);
        self.mesh_vertices = np.vstack((self.mesh_vertices,
            convex_hull.points));
        for face in convex_hull.simplices:
            if is_edge_loop_face(face):
                continue;
            face = self._correct_orientation(joint_center, convex_hull.points, face);
            self.mesh_faces = np.vstack(
                    (self.mesh_faces, face + num_vts));

    def _clip_exterior(self):
        bbox_min, bbox_max = self.original_wire_network.bbox;
        self.mesh_vertices = np.clip(self.mesh_vertices, bbox_min, bbox_max);
        self._clean_up();

    def _enforce_periodic_connectivity(self):
        eps = 1e-6;
        bbox_min, bbox_max = self.original_wire_network.bbox;
        bbox_size = bbox_max - bbox_min;

        vtx_indices = np.arange(len(self.mesh_vertices), dtype=int);

        min_x_vtx = vtx_indices[self.mesh_vertices[:,0] < bbox_min[0] + eps];
        min_y_vtx = vtx_indices[self.mesh_vertices[:,1] < bbox_min[1] + eps];
        min_z_vtx = vtx_indices[self.mesh_vertices[:,2] < bbox_min[2] + eps];
        max_x_vtx = vtx_indices[self.mesh_vertices[:,0] > bbox_max[0] - eps];
        max_y_vtx = vtx_indices[self.mesh_vertices[:,1] > bbox_max[1] - eps];
        max_z_vtx = vtx_indices[self.mesh_vertices[:,2] > bbox_max[2] - eps];

        offset_x = [bbox_size[0], 0.0, 0.0];
        offset_y = [0.0, bbox_size[1], 0.0];
        offset_z = [0.0, 0.0, bbox_size[2]];

        matched_x = self._match_vertices(min_x_vtx, max_x_vtx, offset_x);
        matched_y = self._match_vertices(min_y_vtx, max_y_vtx, offset_y);
        matched_z = self._match_vertices(min_z_vtx, max_z_vtx, offset_z);

        self._correct_inconsistent_faces(*matched_x);
        self._correct_inconsistent_faces(*matched_y);
        self._correct_inconsistent_faces(*matched_z);

    def _match_vertices(self, min_vtx, max_vtx, offset):
        """ Find vertices in max_vtx that matches each vertex in min_vtx.
        """
        assert(len(min_vtx) == len(max_vtx));
        num_pts = len(min_vtx);
        cell_size = np.amax(offset) * 0.001;
        hash_grid = PyMesh.HashGrid.create(cell_size);
        hash_grid.insert_multiple(min_vtx,
                self.mesh_vertices[min_vtx] + offset);

        mapped_min_vtx = np.zeros(num_pts, dtype=int);
        for i,vi in enumerate(max_vtx):
            v = self.mesh_vertices[vi];
            nearby_v_indices = hash_grid.get_items_near_point(v).ravel();
            if len(nearby_v_indices) != 1:
                err_msg = "Inflated mesh cannot be tiled in " +\
                        "[{}] direction.\n".format(offset);
                err_msg += "Found {} vertices to match {}".format(
                        len(nearby_v_indices), v);
                raise RuntimeError(err_msg);
            mapped_min_vtx[i] = nearby_v_indices[0];
        min_vtx = mapped_min_vtx;
        return min_vtx, max_vtx;

    def _correct_inconsistent_faces(self, v_indices_1, v_indices_2):
        """ Given that each item in v_indices_1 matches the corresponding item in
        v_indices_2, make sure faces consisting solely of v_indices_1 match the
        connectivity of faces consisting solely of v_indices_2.
        """
        v_group = np.zeros(len(self.mesh_vertices));
        v_group[v_indices_1] = 1;
        v_group[v_indices_2] = 2;

        v_map_1_to_2 = np.zeros(len(self.mesh_vertices)) - 1;
        v_map_1_to_2[v_indices_1] = v_indices_2;

        face_group = v_group[self.mesh_faces];
        face_group_1 = np.all(face_group == 1, axis=1);
        face_group_2 = np.all(face_group == 2, axis=1);

        faces_in_1 = self.mesh_faces[face_group_1];
        faces_in_1_mapped = v_map_1_to_2[faces_in_1];
        faces_in_1_mapped = faces_in_1_mapped[:, [0, 2, 1]];
        self.mesh_faces[face_group_2] = faces_in_1_mapped;



import os
import sys

import numpy as np
from numpy.linalg import norm
from math import pi, sqrt, radians
from scipy.spatial import ConvexHull

from BoxIntersection import BoxIntersection
from Subdivision import Subdivision
from UniquePointExtractor import UniquePointExtractor
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


    def inflate(self, clean_up=True):
        if not clean_up:
            raise NotImplementedError("Clean up is required for periodic wires");
        super(PeriodicWireInflator, self).inflate(clean_up);
        self._enforce_periodic_connectivity();
        self._clean_up();
        self._subdivide(1);

    def _subdivide(self, num_iterations):
        sub = Subdivision();
        self.mesh_vertices, self.mesh_faces = \
                sub.subdivide(self.mesh_vertices, self.mesh_faces, num_iterations);

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

    def _generate_edge_pipe(self, ei, segment_len, vertex_count):
        edge = self.wire_network.edges[ei];
        if np.any(self.phantom_vertex_map[edge] < 0):
            return [], []
        return super(PeriodicWireInflator, self)._generate_edge_pipe(
                ei, segment_len, vertex_count);

    def _generate_joints(self):
        self._register_edge_loops();
        super(PeriodicWireInflator, self)._generate_joints();

    def _register_edge_loops(self):
        self.grid = PyMesh.HashGrid.create(1e-3);
        for i,edge in self.wire_network.edges:
            end_pts = self.wire_network.vertices[edge];
            loop_0 = self.edge_loops[i, 0, :, :].reshape((-1, 3), order="C");
            loop_1 = self.edge_loops[i, 0, :, :].reshape((-1, 3), order="C");
            self.grid.insert_batch(i, loop_0);
            self.grid.insert_batch(i, loop_1);

    def _look_up_edge_loop_vertex(self, p):
        candidates = self.grid.get_items_near_point(p);
        if len(candidates) == 0:
            return None;
        elif len(candidates) == 1:
            return candidates[0];
        else:
            raise RuntimeError("More than one edge loops occupies this point");

    def _generate_joint(self, idx, num_vts):
        if self.phantom_vertex_map[idx] < 0:
            return 0;

        v = self.wire_network.vertices[idx];
        loop_vertices, phantom_indicator = self._get_incident_edge_loop_vertices(idx);
        loop_vertices = np.vstack([v, loop_vertices]);
        phantom_indicator = [False] + phantom_indicator.tolist();
        hull = ConvexHull(loop_vertices);

        bbox_min, bbox_max = self.original_wire_network.bbox;
        cell_bbox = BoxIntersection(bbox_min, bbox_max);
        intersection_pts = cell_bbox.intersect(hull.points, hull.simplices);

        eps = 1e-8;
        inside = np.logical_and(
                np.all(loop_vertices > bbox_min - eps, axis=1),
                np.all(loop_vertices < bbox_max + eps, axis=1));
        interior_pts = loop_vertices[inside];
        interior_pts = np.clip(interior_pts, bbox_min, bbox_max);

        if len(intersection_pts) > 0:
            vertices = np.vstack((interior_pts, intersection_pts));
        else:
            vertices = interior_pts;
        vertices = UniquePointExtractor.extract(vertices);
        hull = ConvexHull(vertices);
        joint_center = np.mean(vertices, axis=0);

        clipped_loop_vertices = np.clip(loop_vertices, bbox_min, bbox_max);
        index_map = self._map_points(clipped_loop_vertices, hull.points);
        self._register_incident_edge_loop_vertex_indices(idx, index_map, num_vts);

        valid_clipped_loop_vertices = clipped_loop_vertices[
                np.logical_not(phantom_indicator)];
        inverse_map = self._map_points(hull.points, valid_clipped_loop_vertices);

        faces = [];
        for face in hull.simplices:
            face = self._correct_orientation(joint_center, hull.points, face);
            loop_indices = set(range(len(hull.points)));
            for vi in face:
                loop_indices &= set((inverse_map[vi].ravel() - 1)/4);
            if len(loop_indices) == 1:
                continue;
            else:
                faces.append(face);
        faces = np.array(faces, dtype=int);

        self.mesh_vertices.append(hull.points);
        self.mesh_faces.append(faces + num_vts);

        return len(hull.points);

    def _map_points(self, from_pts, to_pts):
        grid = PyMesh.HashGrid.create(1e-8);
        grid.insert_multiple(np.arange(len(to_pts), dtype=int), to_pts);

        indices = [];
        for i,p in enumerate(from_pts):
            candidates = grid.get_items_near_point(p);
            indices.append(candidates);
        return indices;

    def _get_incident_edge_loop_vertices(self, vertex_idx):
        loop_vertices = [];
        phantom_indicator = [];
        for i,e_idx in enumerate(self.wire_network.vertex_edge_neighbors[vertex_idx]):
            edge = self.wire_network.edges[e_idx];
            v_idx = np.arange(2)[edge == vertex_idx];
            loop = self.edge_loops[e_idx, v_idx, :, :].reshape((-1, 3),
                    order="C");
            loop_vertices.append(loop);
            phantom_indicator.append(np.any(self.phantom_vertex_map[edge] < 0));
        loop_vertices = np.vstack(loop_vertices);
        phantom_indicator = np.repeat(phantom_indicator, 4);
        return loop_vertices, phantom_indicator;

    def _register_incident_edge_loop_vertex_indices(self, vertex_idx, index_map,
            num_vts):
        for i, e_idx in enumerate(self.wire_network.vertex_edge_neighbors[vertex_idx]):
            edge = self.wire_network.edges[e_idx];
            if np.any(self.phantom_vertex_map[edge] < 0):
                continue;
            v_idx = np.arange(2)[edge == vertex_idx];
            loop_v_idx = [];
            for vi in range(i*4+1, i*4+4+1):
                assert(len(index_map[vi]) == 1)
                loop_v_idx.append(index_map[vi].ravel()[0]);

            self.edge_loop_indices[e_idx*2+v_idx, :] = \
                    np.array(loop_v_idx, dtype=int) + num_vts;

    def _enforce_periodic_connectivity(self):
        eps = 1e-2;
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

        self.snap_corresponding_vertices(*matched_x);
        self.snap_corresponding_vertices(*matched_y);
        self.snap_corresponding_vertices(*matched_z);

        self._correct_inconsistent_faces(*matched_x);
        self._correct_inconsistent_faces(*matched_y);
        self._correct_inconsistent_faces(*matched_z);

    def _match_vertices(self, min_vtx, max_vtx, offset):
        """ Find vertices in max_vtx that matches each vertex in min_vtx.
        """
        assert(len(min_vtx) == len(max_vtx));
        num_pts = len(min_vtx);
        bbox_min, bbox_max = self.original_wire_network.bbox;
        cell_size = np.amin(bbox_max - bbox_min) * 0.001;
        hash_grid = PyMesh.HashGrid.create(cell_size);
        hash_grid.insert_multiple(min_vtx,
                self.mesh_vertices[min_vtx] + offset);

        mapped_min_vtx = [];
        mapped_max_vtx = [];
        for i,vi in enumerate(max_vtx):
            v = self.mesh_vertices[vi];
            nearby_v_indices = hash_grid.get_items_near_point(v).ravel();
            if len(nearby_v_indices) != 1:
                pass;
                #err_msg = "Inflated mesh cannot be tiled in " +\
                #        "{} direction.\n".format(offset);
                #err_msg += "Found {} vertices to match {}".format(
                #        len(nearby_v_indices), v);
                #raise RuntimeError(err_msg);
            else:
                mapped_min_vtx.append(nearby_v_indices[0]);
                mapped_max_vtx.append(vi);

        min_vtx = np.array(mapped_min_vtx, dtype=int);
        max_vtx = np.array(mapped_max_vtx, dtype=int);
        return min_vtx, max_vtx;

    def snap_corresponding_vertices(self, min_vertices, max_vertices):
        bbox_min, bbox_max = self.original_wire_network.bbox;
        for i in range(len(min_vertices)):
            vi = min_vertices[i];
            vj = max_vertices[i];
            v_min = self.mesh_vertices[vi];
            v_max = self.mesh_vertices[vj];
            direction = np.argmax(v_max - v_min);
            v_mid = 0.5 * (v_min + v_max);
            v_min = np.copy(v_mid);
            v_max = np.copy(v_mid);
            v_min[direction] = bbox_min[direction];
            v_max[direction] = bbox_max[direction];
            self.mesh_vertices[vi] = v_min;
            self.mesh_vertices[vj] = v_max;

    def _correct_inconsistent_faces(self, v_indices_1, v_indices_2):
        """ Given that each item in v_indices_1 matches the corresponding item in
        v_indices_2, make sure faces consisting solely of v_indices_1 match the
        connectivity of faces consisting solely of v_indices_2.
        """
        v_group = np.zeros(len(self.mesh_vertices));
        v_group[v_indices_1] = 1;
        v_group[v_indices_2] = 2;

        v_map_1_to_2 = np.zeros(len(self.mesh_vertices), dtype=int) - 1;
        v_map_1_to_2[v_indices_1] = v_indices_2;

        face_group = v_group[self.mesh_faces];
        face_group_1 = np.all(face_group == 1, axis=1);
        face_group_2 = np.all(face_group == 2, axis=1);
        assert(np.sum(face_group_1) == np.sum(face_group_2));

        faces_in_1 = self.mesh_faces[face_group_1];
        faces_in_1_mapped = v_map_1_to_2[faces_in_1];
        faces_in_1_mapped = faces_in_1_mapped[:, [0, 2, 1]];
        self.mesh_faces[face_group_2] = faces_in_1_mapped;


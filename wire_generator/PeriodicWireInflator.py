import os
import sys

import numpy as np
from numpy.linalg import norm
from math import pi, sqrt, radians
from scipy.spatial import ConvexHull

from BoxIntersection import BoxIntersection
from UniquePointExtractor import UniquePointExtractor
from WireNetwork import WireNetwork
from WireInflator import WireInflator
from WirePattern import WirePattern
from timethis import timethis
from WireWriter import WireWriter

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

        self.__compute_original_bbox();

    def inflate(self, clean_up=True):
        if not clean_up:
            raise NotImplementedError("Clean up is required for periodic wires");
        super(PeriodicWireInflator, self).inflate(clean_up);
        self._enforce_periodic_connectivity();
        self._subdivide(2);

    def __generate_phantom_wire_network(self):
        wire_pattern = WirePattern();
        wire_pattern.set_single_cell_from_wire_network(self.wire_network);
        wire_pattern.tile([3, 3, 3]);
        phantom_wire_network = wire_pattern.wire_network;

        phantom_wire_network.translate(
                self.wire_network.bbox_center - phantom_wire_network.bbox_center);
        self.phantom_wire_network = phantom_wire_network;

    def __generate_phantom_vertex_map(self):
        bbox_min, bbox_max = self.phantom_wire_network.bbox;
        cell_size = np.amin(bbox_max - bbox_min) * 0.001;
        if cell_size == 0:
            cell_size = 1e-3;
        hash_grid = PyMesh.HashGrid.create(cell_size, self.wire_network.dim);
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

    def __compute_original_bbox(self):
        eps = 1e-8;
        bbox_min, bbox_max = self.original_wire_network.bbox;
        self.degenerated_dim = (bbox_max - bbox_min) < eps;
        if np.any(self.degenerated_dim):
            # Input wire network is flat.
            large_value = np.amax(bbox_max - bbox_min);
            bbox_min[self.degenerated_dim] -= large_value;
            bbox_max[self.degenerated_dim] += large_value;
        self.original_bbox_min = bbox_min;
        self.original_bbox_max = bbox_max;

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
        self.grid = PyMesh.HashGrid.create(1e-3, self.wire_network.dim);
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

        eps = 1e-8;
        v = self.wire_network.vertices[idx];
        loop_vertices, phantom_indicator = self._get_incident_edge_loop_vertices(idx);
        loop_vertices = np.vstack([v, loop_vertices]);
        phantom_indicator = [False] + phantom_indicator.tolist();
        hull = ConvexHull(loop_vertices);

        bbox_min = self.original_bbox_min;
        bbox_max = self.original_bbox_max;
        cell_bbox = BoxIntersection(bbox_min, bbox_max);
        intersection_pts = cell_bbox.intersect(hull.points, hull.simplices);

        inside = np.logical_and(
                np.all(loop_vertices > bbox_min - eps, axis=1),
                np.all(loop_vertices < bbox_max + eps, axis=1));
        interior_pts = loop_vertices[inside];
        interior_pts = np.clip(interior_pts, bbox_min, bbox_max);

        if len(intersection_pts) > 0:
            vertices = np.vstack((interior_pts, intersection_pts));
        else:
            vertices = interior_pts;
        assert(len(vertices) > 0);
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
        grid = PyMesh.HashGrid.create(1e-8, self.wire_network.dim);
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


    def enforce_single_axis_periodicity(self, axis):
        if self.degenerated_dim[axis]: return;

        tol = 1e-12;
        bbox_min = self.original_bbox_min;
        bbox_max = self.original_bbox_max;

        vertices = self.mesh_vertices;
        faces = self.mesh_faces;
        num_vertices = len(self.mesh_vertices);

        v_on_min_boundary = vertices[:, axis] <= bbox_min[axis] + tol;
        v_on_max_boundary = vertices[:, axis] >= bbox_max[axis] - tol;
        assert(np.any(v_on_min_boundary));
        assert(np.any(v_on_max_boundary));

        offset = np.zeros(self.wire_network.dim);
        offset[axis] = bbox_max[axis] - bbox_min[axis];

        f_on_min_boundary = np.all(v_on_min_boundary[faces], axis=1);
        f_on_max_boundary = np.all(v_on_max_boundary[faces], axis=1);
        assert(np.any(f_on_min_boundary));
        assert(np.any(f_on_max_boundary));

        added_vertices = vertices[v_on_min_boundary] + offset;
        from_index = np.arange(num_vertices, dtype=int)[v_on_min_boundary];
        to_index = np.arange(len(added_vertices), dtype=int) + num_vertices;
        vertex_map = {i:j for i,j in zip(from_index, to_index)};
        index_lookup = lambda i: vertex_map[i];
        added_faces = [map(index_lookup, face) for face in faces[f_on_min_boundary]];
        added_faces = np.array(added_faces)[:,[1,0,2]];

        vertices = np.vstack((vertices, added_vertices));
        faces = faces[np.logical_not(f_on_max_boundary)];
        faces = np.vstack((faces, added_faces));

        self.mesh_vertices = vertices;
        self.mesh_faces = faces;

    def _enforce_periodic_connectivity(self):
        self.enforce_single_axis_periodicity(0);
        self._clean_up();
        self.enforce_single_axis_periodicity(1);
        self._clean_up();
        self.enforce_single_axis_periodicity(2);
        self._clean_up();

    def write_debug_wires(self, filename, vertices, edges):
        writer = WireWriter(filename);
        writer.write(vertices, edges);


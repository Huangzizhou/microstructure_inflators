import os
import sys

import numpy as np
from numpy.linalg import norm
from math import pi, sqrt
from scipy.spatial import ConvexHull

from WireNetwork import WireNetwork
from Subdivision import Subdivision
from timethis import timethis

import PyMeshSetting
import PyMesh
import PyMeshUtils

class WireInflator(object):
    def __init__(self, wire_network):
        self.wire_network = wire_network;
        if self.wire_network.dim == 2:
            raise NotImplementedError(
                    "Inflating 2D wire network is not supported yet.");

    @timethis
    def inflate(self, clean_up=True, subdivide_order=1):
        self._compute_thickness();
        self._compute_min_edge_angles();
        self._compute_edge_end_loops();
        self._generate_joints();
        self._generate_edge_pipes();
        if clean_up:
            self._clean_up();
        self._subdivide(subdivide_order);

    @property
    def mesh(self):
        factory = PyMesh.MeshFactory();
        factory.load_data(
                self.mesh_vertices.ravel(order="C"),
                self.mesh_faces.ravel(order="C"),
                np.zeros(0),
                3, 3, 4);
        mesh = factory.create_shared();
        mesh.add_attribute("source_wire_id");
        mesh.set_attribute("source_wire_id", self.source_wire_id);

        if "vertex_orbit_index" in self.wire_network.attributes:
            indices = self.wire_network.attributes["vertex_orbit_index"].ravel();
            source_id_mask = self.source_wire_id > 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("vertex_orbit");
            mesh.set_attribute("vertex_orbit", source_index);

        if "edge_orbit_index" in self.wire_network.attributes:
            indices = self.wire_network.attributes["edge_orbit_index"].ravel();
            source_id_mask = self.source_wire_id < 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[-self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("edge_orbit");
            mesh.set_attribute("edge_orbit", source_index);

        return mesh;

    @timethis
    def save(self, mesh_file):
        mesh = self.mesh;
        writer = PyMesh.MeshWriter.create_writer(mesh_file);
        writer.with_attribute("source_wire_id");
        if mesh.has_attribute("vertex_orbit"):
            writer.with_attribute("vertex_orbit");
        if mesh.has_attribute("edge_orbit"):
            writer.with_attribute("edge_orbit");

        writer.write_mesh(mesh);

    @timethis
    def _compute_thickness(self):
        """ Compute self.thickness parameter.  It is a num_edges x 2 matrix,
        where each row corresponds to the thickness at each end of the edge.
        """
        self.thickness = np.zeros((self.wire_network.num_edges, 2));
        if "vertex_thickness" in self.wire_network.attributes:
            self.__init_vertex_thickness(self.wire_network.attributes["vertex_thickness"]);
        elif "edge_thickness" in self.wire_network.attributes:
            self.__init_edge_thickness(self.wire_network.attributes["edge_thickness"]);
        else:
            raise RuntimeError("Wire thickness not specified through attributes");

        min_thickness = np.amin(self.thickness);
        if min_thickness < 1e-3:
            raise RuntimeError("Extremely thin thickness detected: {}".format(
                min_thickness));

    def __init_vertex_thickness(self, v_thickness):
        for i,e in enumerate(self.wire_network.edges):
            self.thickness[i, :] = v_thickness[e].ravel();

    def __init_edge_thickness(self, e_thickness):
        assert(len(e_thickness) == len(self.thickness));
        for i,e in enumerate(self.wire_network.edges):
            self.thickness[i, :] = e_thickness[i];

    @timethis
    def _compute_min_edge_angles(self):
        self.min_angles = np.zeros(self.wire_network.num_vertices);
        for i in range(self.wire_network.num_vertices):
            angle = self._compute_min_edge_angle(i);
            self.min_angles[i] = angle;

    @timethis
    def _compute_edge_end_loops(self):
        vertex_offsets = self._compute_vertex_offsets();
        self.edge_loops = np.zeros((self.wire_network.num_edges, 2, 4, 3));
        for i,edge in enumerate(self.wire_network.edges):
            angles = self.min_angles[edge];
            thickness = self.thickness[i];
            #offsets = 0.5 * sqrt(2) * thickness / np.tan(angles / 2.0) +1e-6;
            offsets = vertex_offsets[edge];
            edge_dir, perp1_dir, perp2_dir = self._generate_frame(edge);
            v0 = self.wire_network.vertices[edge[0]];
            v1 = self.wire_network.vertices[edge[1]];

            loop_0 = np.array([
                    v0+offsets[0]*edge_dir+0.5*thickness[0]*(-perp1_dir-perp2_dir),
                    v0+offsets[0]*edge_dir+0.5*thickness[0]*(-perp1_dir+perp2_dir),
                    v0+offsets[0]*edge_dir+0.5*thickness[0]*( perp1_dir+perp2_dir),
                    v0+offsets[0]*edge_dir+0.5*thickness[0]*( perp1_dir-perp2_dir) ]);
            loop_1 = np.array([
                    v1-offsets[1]*edge_dir+0.5*thickness[1]*(-perp1_dir-perp2_dir),
                    v1-offsets[1]*edge_dir+0.5*thickness[1]*(-perp1_dir+perp2_dir),
                    v1-offsets[1]*edge_dir+0.5*thickness[1]*( perp1_dir+perp2_dir),
                    v1-offsets[1]*edge_dir+0.5*thickness[1]*( perp1_dir-perp2_dir) ]);

            self.edge_loops[i, 0, :, :] = loop_0;
            self.edge_loops[i, 1, :, :] = loop_1;

    def _compute_vertex_offsets(self):
        vertex_thickness = np.zeros(self.wire_network.num_vertices);
        for i, edge in enumerate(self.wire_network.edges):
            thickness = self.thickness[i];
            vertex_thickness[edge[0]] = max(vertex_thickness[edge[0]], thickness[0]);
            vertex_thickness[edge[1]] = max(vertex_thickness[edge[1]], thickness[1]);

        offsets = np.zeros(self.wire_network.num_vertices);
        offsets = 0.5 * sqrt(2) * vertex_thickness /\
                np.tan(self.min_angles / 2.0) + 1e-6;

        return offsets;


    @timethis
    def _generate_joints(self):
        self.mesh_vertices = [];
        self.mesh_faces = [];
        self.edge_loop_indices = np.zeros((self.edge_loops.shape[0]*2, 4),
                dtype=int);

        # Source wire_id assigns an index per face.
        # positive indices is vertex index.  Negative indices is edge index.
        # 0 should not be used.
        self.source_wire_id = [];

        num_vts = 0;
        num_faces = 0;
        for i in range(self.wire_network.num_vertices):
            num_added_vts = self._generate_joint(i, num_vts);
            num_added_faces = len(self.mesh_faces) - num_faces;
            num_vts += num_added_vts;
            num_faces += num_added_faces;
            self.source_wire_id.append([i+1]*num_added_faces);

        self.mesh_vertices = np.vstack(self.mesh_vertices);
        self.mesh_faces = np.vstack(self.mesh_faces);
        self.source_wire_id = np.hstack(self.source_wire_id);
        assert(len(self.source_wire_id) == len(self.mesh_faces));

    @timethis
    def _generate_edge_pipes(self):
        extra_vertices = [];
        extra_faces = [];
        extra_wire_id = [];
        vertex_count = len(self.mesh_vertices);

        for i,edge in enumerate(self.wire_network.edges):
            segment_len = 2 * np.min(self.thickness[i]);
            new_v, new_f = self._generate_edge_pipe(i, segment_len,
                    vertex_count);

            extra_vertices += new_v;
            extra_faces += new_f;
            extra_wire_id += [-i-1] * len(new_f);
            vertex_count += len(new_v) * 4;

        self.mesh_vertices = np.vstack([self.mesh_vertices] + extra_vertices);
        self.mesh_faces = np.vstack([self.mesh_faces] + extra_faces);
        self.source_wire_id = np.hstack([self.source_wire_id] + extra_wire_id);
        assert(len(self.mesh_vertices) == vertex_count);
        assert(len(self.source_wire_id) == len(self.mesh_faces));

    @timethis
    def _generate_edge_pipe(self, ei, segment_len, vertex_count):
        edge = self.wire_network.edges[ei];
        v0 = self.wire_network.vertices[edge[0]];
        v1 = self.wire_network.vertices[edge[1]];
        edge_len = norm(v1 - v0);
        num_segments = int(np.round(edge_len / segment_len));
        num_segments = max(1, num_segments);

        loop_1_idx = self.edge_loop_indices[2*ei  , :];
        loop_2_idx = self.edge_loop_indices[2*ei+1, :];
        loop_1_vts = self.mesh_vertices[loop_1_idx];
        loop_2_vts = self.mesh_vertices[loop_2_idx];

        loop_idx = [loop_1_idx];
        extra_vertices = [];
        for i in range(1, num_segments):
            frac = float(i) / float(num_segments);
            loop_vtx = loop_1_vts * (1.0 - frac) + loop_2_vts * frac;
            loop_idx.append(range(vertex_count+(i-1)*4, vertex_count+i*4));
            extra_vertices.append(loop_vtx);
        loop_idx.append(loop_2_idx);

        faces = [];
        for i in range(num_segments):
            l1_idx = loop_idx[i];
            l2_idx = loop_idx[i+1];
            faces.append([l1_idx[0], l1_idx[1], l2_idx[0]]);
            faces.append([l2_idx[0], l1_idx[1], l2_idx[1]]);
            faces.append([l2_idx[1], l1_idx[1], l1_idx[2]]);
            faces.append([l2_idx[1], l1_idx[2], l2_idx[2]]);
            faces.append([l2_idx[2], l1_idx[2], l1_idx[3]]);
            faces.append([l2_idx[2], l1_idx[3], l2_idx[3]]);
            faces.append([l2_idx[3], l1_idx[3], l1_idx[0]]);
            faces.append([l2_idx[3], l1_idx[0], l2_idx[0]]);

        return extra_vertices, faces;

    @timethis
    def _compute_min_edge_angle(self, idx):
        min_angle = pi;
        v = self.wire_network.vertices[idx];
        neighbors = self.wire_network.vertex_neighbors[idx];
        num_neighbors = len(neighbors);
        for i in range(num_neighbors):
            vi = self.wire_network.vertices[neighbors[i]];
            ei = vi - v;
            for j in range(i+1, num_neighbors):
                vj = self.wire_network.vertices[neighbors[j]];
                ej = vj - v;
                angle = np.arctan2(norm(np.cross(ei, ej)), np.dot(ei, ej));
                if angle < 0.0:
                    angle += pi;
                min_angle = min(min_angle, angle);
        return min_angle;

    @timethis
    def _generate_frame(self, edge):
        edge_dir = self.wire_network.vertices[edge[1]] -\
                self.wire_network.vertices[edge[0]];
        edge_dir /= norm(edge_dir);
        z_dir = [0.0, 0.0, 1.0];
        offset_dir_1 = np.cross(edge_dir, z_dir);
        offset_dir_1_len = norm(offset_dir_1);
        if offset_dir_1_len < 1e-6:
            offset_dir_1 = np.array([0.0, 1.0, 0.0]);
        else:
            offset_dir_1 /= offset_dir_1_len;
        offset_dir_2 = np.cross(offset_dir_1, edge_dir);
        return edge_dir, offset_dir_1, offset_dir_2;

    @timethis
    def _generate_joint(self, idx, num_vts):
        v = self.wire_network.vertices[idx];
        loop_vertices = [v];
        for i,e_idx in enumerate(self.wire_network.vertex_edge_neighbors[idx]):
            edge = self.wire_network.edges[e_idx];
            if edge[0] == idx:
                v_idx = 0;
            else:
                v_idx = 1;
            loop = self.edge_loops[e_idx, v_idx, :, :].reshape((-1,3), order="C");
            loop_vertices.append(loop);
            self.edge_loop_indices[e_idx * 2 + v_idx, :] = \
                    np.arange(i*4, i*4+4) + num_vts + 1;

        loop_vertices = np.vstack(loop_vertices);
        joint_center = np.mean(loop_vertices, axis=0);
        convex_hull = ConvexHull(loop_vertices);
        self.mesh_vertices.append(convex_hull.points);

        for face in convex_hull.simplices:
            # The following check if a face is made entirely of edge loop
            # vertices.  If this is the case, it should be removed since edge
            # loops would be connected by pipes.
            loop_idx = [(vi-1) / 4 for vi in face];
            if np.max(loop_idx) == np.min(loop_idx):
                continue;
            face = self._correct_orientation(joint_center, convex_hull.points, face);
            self.mesh_faces.append(face + num_vts);

        return len(convex_hull.points);

    @timethis
    def _correct_orientation(self, center, points, face):
        vts = np.array(points)[face];
        face_center = np.mean(vts, axis=0);
        out_dir = face_center - center;
        normal = np.cross(vts[1] - vts[0], vts[2] - vts[0]);
        sign = np.dot(out_dir, normal);
        if sign < 0:
            return face[[0, 2, 1]];
        else:
            return face;

    @timethis
    def _clean_up(self, track_source_id = True):
        ## Remove obtuse triangles
        #obtuse_triangle_removal = PyMeshUtils.ObtuseTriangleRemoval(
        #        self.mesh_vertices, self.mesh_faces);
        #obtuse_triangle_removal.run(179.0);
        #self.mesh_vertices = obtuse_triangle_removal.get_vertices();
        #self.mesh_faces = obtuse_triangle_removal.get_faces();

        # Remove duplicated vertices
        duplicated_vertex_removal = PyMeshUtils.DuplicatedVertexRemoval(
                self.mesh_vertices, self.mesh_faces);
        duplicated_vertex_removal.run(1e-3);
        self.mesh_vertices = duplicated_vertex_removal.get_vertices();
        self.mesh_faces = duplicated_vertex_removal.get_faces();

        # Collapse short edges
        edge_remover = PyMeshUtils.ShortEdgeRemoval(self.mesh_vertices,
                self.mesh_faces);
        edge_remover.run(np.amin(1e-12, 0.1 * np.amin(self.thickness)));
        self.mesh_vertices = edge_remover.get_vertices();
        self.mesh_faces = edge_remover.get_faces();
        if track_source_id:
            self.face_indices = edge_remover.get_face_indices().ravel();
            self.source_wire_id = self.source_wire_id[self.face_indices];

        # Remove isolated vertices
        unique_indices = np.unique(self.mesh_faces.ravel());
        index_map = np.zeros(len(self.mesh_vertices), dtype=int);
        index_map[unique_indices] = np.arange(len(unique_indices), dtype=int);
        self.mesh_vertices = self.mesh_vertices[unique_indices];
        self.mesh_faces = index_map[self.mesh_faces];

        # Edge loop indices is no longer valid.
        self.edge_loop_indices =None;

    @timethis
    def _subdivide(self, num_iterations):
        sub = Subdivision();
        self.mesh_vertices, self.mesh_faces, face_indices= \
                sub.subdivide(self.mesh_vertices, self.mesh_faces, num_iterations);
        self.source_wire_id = self.source_wire_id[face_indices];


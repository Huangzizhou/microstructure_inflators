import os
import sys

import numpy as np
from numpy.linalg import norm
from math import pi
from scipy.spatial import ConvexHull

from WireNetwork import WireNetwork
from timethis import timethis

# Update path to import PyMesh
py_mesh_path = os.environ.get("PYMESH_PATH");
if py_mesh_path == None:
    raise ImportError("Please set PYMESH_PATH to the correct lib path.");
sys.path.append(os.path.join(py_mesh_path, "lib"));
sys.path.append(os.path.join(py_mesh_path, "swig"));
import PyMesh

class WireInflator(object):
    def __init__(self, wire_network):
        self.wire_network = wire_network;
        if self.wire_network.dim == 2:
            NotImplementedError("2D wire inflation is not supported.");

    @timethis
    def inflate(self, thickness):
        self.thickness = thickness;
        self.__compute_min_edge_angles();
        self.__compute_edge_end_loops();
        self.__generate_joints();
        self.__generate_edge_pipes();

    @timethis
    def save(self, mesh_file):
        writer = PyMesh.MeshWriter.create_writer(mesh_file);
        writer.write(
                self.mesh_vertices.ravel(order="C"),
                self.mesh_faces.ravel(order="C"),
                np.zeros(0),
                3, 3, 4);

    @timethis
    def __compute_min_edge_angles(self):
        self.min_angles = np.zeros(self.wire_network.num_vertices);
        for i in range(self.wire_network.num_vertices):
            angle = self.__compute_min_edge_angle(i);
            self.min_angles[i] = angle;

    @timethis
    def __compute_edge_end_loops(self):
        self.edge_loops = np.zeros((self.wire_network.num_edges, 2, 4, 3));
        for i,edge in enumerate(self.wire_network.edges):
            angles = self.min_angles[edge];
            offsets = 0.5 * self.thickness / np.tan(angles / 2.0) + self.thickness;
            edge_dir, perp1_dir, perp2_dir = self.__generate_frame(edge);
            v0 = self.wire_network.vertices[edge[0]];
            v1 = self.wire_network.vertices[edge[1]];

            loop_0 = [
                    v0+offsets[0]*edge_dir+0.5*self.thickness*(-perp1_dir-perp2_dir),
                    v0+offsets[0]*edge_dir+0.5*self.thickness*(-perp1_dir+perp2_dir),
                    v0+offsets[0]*edge_dir+0.5*self.thickness*( perp1_dir+perp2_dir),
                    v0+offsets[0]*edge_dir+0.5*self.thickness*( perp1_dir-perp2_dir) ]
            loop_1 = [
                    v1-offsets[1]*edge_dir+0.5*self.thickness*(-perp1_dir-perp2_dir),
                    v1-offsets[1]*edge_dir+0.5*self.thickness*(-perp1_dir+perp2_dir),
                    v1-offsets[1]*edge_dir+0.5*self.thickness*( perp1_dir+perp2_dir),
                    v1-offsets[1]*edge_dir+0.5*self.thickness*( perp1_dir-perp2_dir) ]
            self.edge_loops[i, 0, :, :] = loop_0;
            self.edge_loops[i, 1, :, :] = loop_1;

    @timethis
    def __generate_joints(self):
        self.mesh_vertices = np.zeros((0, 3), dtype=float);
        self.mesh_faces = np.zeros((0, 3), dtype=int);
        self.edge_loop_indices = np.zeros((self.edge_loops.shape[0]*2, 4),
                dtype=int);
        for i in range(self.wire_network.num_vertices):
            self.__generate_joint(i);

    @timethis
    def __generate_edge_pipes(self):
        for i,edge in enumerate(self.wire_network.edges):
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

    @timethis
    def __compute_min_edge_angle(self, idx):
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
    def __generate_frame(self, edge):
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
    def __generate_joint(self, idx):
        num_vts = len(self.mesh_vertices);
        v = self.wire_network.vertices[idx];
        loop_vertices = v.reshape((-1,3));
        for i,e_idx in enumerate(self.wire_network.vertex_edge_neighbors[idx]):
            edge = self.wire_network.edges[e_idx];
            v_idx = np.where(np.array(edge).ravel() == idx)[0][0];
            loop = self.edge_loops[e_idx, v_idx, :, :].reshape((-1,3), order="C");
            loop_vertices = np.vstack((loop_vertices, loop));
            self.edge_loop_indices[e_idx * 2 + v_idx, :] = \
                    np.arange(i*4, i*4+4) + num_vts + 1;

        joint_center = np.mean(loop_vertices, axis=0);
        convex_hull = ConvexHull(loop_vertices);
        self.mesh_vertices = np.vstack((self.mesh_vertices,
            convex_hull.points));
        for face in convex_hull.simplices:
            loop_idx = [(vi-1) / 4 for vi in face];
            if np.max(loop_idx) == np.min(loop_idx):
                continue;
            face = self.__correct_orientation(joint_center, convex_hull.points, face);
            self.mesh_faces = np.vstack(
                    (self.mesh_faces, face + num_vts));

    def __correct_orientation(self, center, points, face):
        vts = np.array(points)[face];
        face_center = np.mean(vts, axis=0);
        out_dir = face_center - center;
        normal = np.cross(vts[1] - vts[0], vts[2] - vts[0]);
        sign = np.dot(out_dir, normal);
        if sign < 0:
            return face[[0, 2, 1]];
        else:
            return face;



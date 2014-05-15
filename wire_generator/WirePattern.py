import numpy as np
import os
import sys

import PyMeshSetting
import PyMesh

from WireNetwork import WireNetwork
from timethis import timethis

class WirePattern(object):
    def __init__(self):
        self.x_tile_dir = np.array([1.0, 0.0, 0.0]);
        self.y_tile_dir = np.array([0.0, 1.0, 0.0]);
        self.z_tile_dir = np.array([0.0, 0.0, 1.0]);

    def set_single_cell(self, vertices, edges):
        bbox_min = np.amin(vertices, axis=0);
        bbox_max = np.amax(vertices, axis=0);
        centroid = 0.5 * (bbox_min + bbox_max);
        self.pattern_vertices = np.array(vertices) - centroid;
        self.pattern_edges = np.array(edges, dtype=int);
        self.pattern_bbox_min = bbox_min - centroid;
        self.pattern_bbox_max = bbox_max - centroid;
        self.pattern_bbox_size = self.pattern_bbox_max - self.pattern_bbox_min;

    def set_single_cell_from_wire_network(self, network):
        self.set_single_cell(network.vertices, network.edges);
        self.pattern_attributes = network.attributes;

    @timethis
    def tile(self, reps):
        self.wire_vertices = np.zeros((0, 3));
        self.wire_edges = np.zeros((0, 2), dtype=int);
        self.pattern_vertex_map = [];
        self.pattern_edge_map = [];
        for i in range(reps[0]):
            x_inc = self.x_tile_dir * self.pattern_bbox_size[0] * i;
            for j in range(reps[1]):
                y_inc = self.y_tile_dir * self.pattern_bbox_size[1] * j;
                for k in range(reps[2]):
                    z_inc = self.z_tile_dir * self.pattern_bbox_size[2] * k;

                    base_idx = self.wire_vertices.shape[0];
                    vertices = self.pattern_vertices +\
                            x_inc + y_inc + z_inc;
                    self.wire_vertices = np.vstack(
                            (self.wire_vertices, vertices));
                    self.pattern_vertex_map += range(len(self.pattern_vertices));
                    self.pattern_edge_map += range(len(self.pattern_edges));

                    edges = self.pattern_edges + base_idx;
                    self.wire_edges = np.vstack(
                            (self.wire_edges, edges));
        self.__remove_duplicated_vertices();
        self.__remove_duplicated_edges();
        self.__center_at_origin();

    @timethis
    def tile_hex_mesh(self, mesh):
        dim = mesh.get_dim();
        if not mesh.has_attribute("voxel_centroid"):
            mesh.add_attribute("voxel_centroid");
        centroids = mesh.get_attribute("voxel_centroid");
        centroids = centroids.reshape((-1, dim), order="C");
        mesh_vertices = mesh.get_vertices().reshape((-1, dim), order="C");

        wire_vertices = [];
        wire_edges = [];
        num_pattern_vertices = len(self.pattern_vertices);
        for i in range(mesh.get_num_voxels()):
            voxel = mesh.get_voxel(i).ravel();
            voxel_vertices = mesh_vertices[voxel];
            voxel_bbox_min = np.amin(voxel_vertices, axis=0);
            voxel_bbox_max = np.amax(voxel_vertices, axis=0);
            scale = np.divide(voxel_bbox_max - voxel_bbox_min,
                    self.pattern_bbox_size);

            base_idx = num_pattern_vertices * i;
            c = centroids[i];

            vertices = self.pattern_vertices * scale + c;
            edges = self.pattern_edges + base_idx;

            wire_vertices.append(vertices);
            wire_edges.append(edges);

        self.wire_vertices = np.vstack(wire_vertices);
        self.wire_edges = np.vstack(wire_edges);
        self.__remove_duplicated_vertices();
        self.__remove_duplicated_edges();
        self.__center_at_origin();

    @timethis
    def __remove_duplicated_vertices(self):
        cell_size = np.amin(self.pattern_bbox_size) * 0.01;
        hash_grid = PyMesh.HashGrid.create(cell_size);
        hash_grid.insert_multiple(
                np.arange(len(self.wire_vertices), dtype=int),
                self.wire_vertices);

        vertices = [];
        v_map = np.arange(len(self.wire_vertices), dtype=int);
        for i,v in enumerate(self.wire_vertices):
            if v_map[i] != i: continue;
            nearby_v_indices = hash_grid.get_items_near_point(v);
            v_map[nearby_v_indices] = len(vertices);
            vertices.append(v);

        self.wire_vertices = np.vstack(vertices);
        self.wire_edges = v_map[self.wire_edges];

        pattern_v_map = np.zeros(len(self.wire_vertices), dtype=int);
        for i,vi in enumerate(v_map):
            pattern_v_map[vi] = self.pattern_vertex_map[i];
        self.pattern_vertex_map = pattern_v_map;

    @timethis
    def __remove_duplicated_edges(self):
        edges = [tuple(sorted(edge)) for edge in self.wire_edges];
        self.wire_edges = np.unique(edges);
        edge_map = {tuple(edge):i for i,edge in enumerate(self.wire_edges)};
        edge_index_map = [edge_map[edge] for edge in edges];

        pattern_e_map = np.zeros(len(self.wire_edges), dtype=int);
        for i,ei in enumerate(edge_index_map):
            pattern_e_map[ei] = self.pattern_edge_map[i];
        self.pattern_edge_map = pattern_e_map;

    @timethis
    def __center_at_origin(self):
        bbox_min = np.amin(self.wire_vertices, axis=0);
        bbox_max = np.amax(self.wire_vertices, axis=0);
        bbox_center = 0.5 * (bbox_min + bbox_max);
        self.wire_vertices = self.wire_vertices - bbox_center;

    @timethis
    def __copy_attributes(self, wire_network):
        """ Copy any vertex attributes from single cell to tiled wire network.
        """
        if hasattr(self, "pattern_attributes"):
            for attr_name in self.pattern_attributes:
                value = self.pattern_attributes[attr_name];
                if len(value) == len(self.pattern_vertices):
                    tiled_value = value[self.pattern_vertex_map];
                    assert(len(tiled_value) == len(self.wire_vertices));
                    wire_network.attributes.add(attr_name, tiled_value);
                elif len(value) == len(self.pattern_edges):
                    tiled_value = value[self.pattern_edge_map];
                    assert(len(tiled_value) == len(self.wire_edges));
                    wire_network.attributes.add(attr_name, tiled_value);

    @property
    @timethis
    def wire_network(self):
        network = WireNetwork();
        network.load(self.wire_vertices, self.wire_edges);
        self.__copy_attributes(network);
        return network;


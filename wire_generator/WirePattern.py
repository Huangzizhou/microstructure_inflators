import numpy as np
import copy
import os
import sys

import PyMeshSetting
import PyMesh
import PyMeshUtils

from WireNetwork import WireNetwork
from timethis import timethis

class WirePattern(object):
    def __init__(self):
        self.x_tile_dir = np.array([1.0, 0.0, 0.0]);
        self.y_tile_dir = np.array([0.0, 1.0, 0.0]);
        self.z_tile_dir = np.array([0.0, 0.0, 1.0]);

    def set_single_cell(self, vertices, edges):
        wire_network = WireNetwork();
        wire_network.load(vertices, edges);
        self.set_single_cell_from_wire_network(wire_network);

    def set_single_cell_from_wire_network(self, network):
        network.offset(-network.bbox_center);
        self.pattern = network;
        self.pattern_bbox_min = np.amin(self.pattern.vertices, axis=0);
        self.pattern_bbox_max = np.amax(self.pattern.vertices, axis=0);
        self.pattern_bbox_size = self.pattern_bbox_max - self.pattern_bbox_min;

    def tile_once(self, base_idx, wire_network, modifiers, **kwargs):
        for modifier in modifiers:
            modifier.modify(wire_network, **kwargs);

        self.wire_vertices.append(wire_network.vertices);
        self.wire_edges.append(wire_network.edges + base_idx);
        for key,val in wire_network.attributes.iteritems():
            self.wire_attributes[key] = self.wire_attributes.get(key, []) + list(val);

    @timethis
    def tile(self, reps, modifiers=[]):
        self.wire_vertices = [];
        self.wire_edges = [];
        self.wire_attributes = {};

        base_idx = 0;
        for i in range(reps[0]):
            x_inc = self.x_tile_dir * self.pattern_bbox_size[0] * i;
            for j in range(reps[1]):
                y_inc = self.y_tile_dir * self.pattern_bbox_size[1] * j;
                for k in range(reps[2]):
                    z_inc = self.z_tile_dir * self.pattern_bbox_size[2] * k;

                    pattern = copy.deepcopy(self.pattern);
                    pattern.vertices += x_inc + y_inc + z_inc;
                    self.tile_once(base_idx, pattern, modifiers);
                    base_idx += pattern.num_vertices;
        self.wire_vertices = np.vstack(self.wire_vertices);
        self.wire_edges = np.vstack(self.wire_edges);
        self.__remove_duplicated_vertices();
        self.__remove_duplicated_edges();
        self.__center_at_origin();
        self.__apply_vertex_offset();

    @timethis
    def tile_hex_mesh(self, mesh, modifiers=[]):
        dim = mesh.get_dim();
        if not mesh.has_attribute("voxel_centroid"):
            mesh.add_attribute("voxel_centroid");
        centroids = mesh.get_attribute("voxel_centroid");
        centroids = centroids.reshape((-1, dim), order="C");
        mesh_vertices = mesh.get_vertices().reshape((-1, dim), order="C");

        mesh_attributes = {};
        for name in mesh.get_attribute_names():
            attr = mesh.get_attribute(name);
            if len(attr) == mesh.get_num_vertices():
                attr = PyMeshUtils.convert_to_voxel_attribute_from_name(name);
            mesh_attributes[name] = attr.ravel();

        self.wire_vertices = [];
        self.wire_edges = [];
        self.wire_attributes = {};
        base_idx = 0;
        for i in range(mesh.get_num_voxels()):
            voxel = mesh.get_voxel(i).ravel();
            voxel_vertices = mesh_vertices[voxel];
            attr_dict = {
                    name:val[i] for name,val in mesh_attributes.iteritems() };

            pattern = copy.deepcopy(self.pattern);
            pattern.vertices = self.trilinear_interpolate(
                    pattern.vertices, voxel_vertices);
            self.tile_once(base_idx, pattern, modifiers, **attr_dict);
            base_idx += pattern.num_vertices;

        self.wire_vertices = np.vstack(self.wire_vertices);
        self.wire_edges = np.vstack(self.wire_edges);
        self.__remove_duplicated_vertices();
        self.__remove_duplicated_edges();
        self.__center_at_origin();
        self.__apply_vertex_offset();

    @timethis
    def trilinear_interpolate(self, vertices, corners):
        """ Trilinear interpolate vertices inside of the hex specified by
        corners.  The order of corners should be consistent with the hex node
        order specified by GMSH:
        http://www.geuz.org/gmsh/doc/texinfo/gmsh.html#Node-ordering
        """
        assert(len(corners) == 8);
        shape_functions = [
                lambda v : (1.0-v[0]) * (1.0-v[1]) * (1.0-v[2]),
                lambda v : (    v[0]) * (1.0-v[1]) * (1.0-v[2]),
                lambda v : (    v[0]) * (    v[1]) * (1.0-v[2]),
                lambda v : (1.0-v[0]) * (    v[1]) * (1.0-v[2]),
                lambda v : (1.0-v[0]) * (1.0-v[1]) * (    v[2]),
                lambda v : (    v[0]) * (1.0-v[1]) * (    v[2]),
                lambda v : (    v[0]) * (    v[1]) * (    v[2]),
                lambda v : (1.0-v[0]) * (    v[1]) * (    v[2]),
                ];
        bbox_min = np.amin(vertices, axis=0);
        bbox_max = np.amax(vertices, axis=0);
        bbox_size = bbox_max - bbox_min;
        vertices = (vertices - bbox_min) / bbox_size;
        trilinear_coords = np.array([
            [f(v) for f in shape_functions]
            for v in vertices ]);
        vertices = np.dot(trilinear_coords, corners);
        return vertices;

    @timethis
    def __remove_duplicated_vertices(self):
        num_vertices = len(self.wire_vertices);
        cell_size = np.amin(self.pattern_bbox_size) * 0.01;
        hash_grid = PyMesh.HashGrid.create(cell_size);
        hash_grid.insert_multiple(
                np.arange(num_vertices, dtype=int),
                self.wire_vertices);

        vertices = [];
        v_map = np.arange(num_vertices, dtype=int);
        for i,v in enumerate(self.wire_vertices):
            if v_map[i] != i: continue;
            nearby_v_indices = hash_grid.get_items_near_point(v);
            v_map[nearby_v_indices] = len(vertices);
            vertices.append(v);

        self.wire_vertices = np.vstack(vertices);
        self.wire_edges = v_map[self.wire_edges];

        # Update vertex attributes
        for key,val in self.wire_attributes.iteritems():
            if len(val) % num_vertices == 0:
                per_vertex_size = len(val) / num_vertices;
                attr_val = [[] for i in range(len(self.wire_vertices))];
                for i,vi in enumerate(v_map):
                    attr_val[vi].append(val[
                        i*per_vertex_size:(i+1)*per_vertex_size]);
                attr_val = [np.mean(val, axis=0) for val in attr_val];
                self.wire_attributes[key] = np.array(attr_val);


                #attr_val = np.zeros(len(self.wire_vertices))
                #for i,vi in enumerate(v_map):
                #    attr_val[vi] = val[i];
                #self.wire_attributes[key] = attr_val;

    @timethis
    def __remove_duplicated_edges(self):
        num_edges = len(self.wire_edges);
        edges = [tuple(sorted(edge)) for edge in self.wire_edges];
        self.wire_edges = np.unique(edges);
        edge_map = {tuple(edge):i for i,edge in enumerate(self.wire_edges)};
        edge_index_map = [edge_map[edge] for edge in edges];

        for key,val in self.wire_attributes.iteritems():
            if len(val) % num_edges == 0:
                per_edge_size = len(val) / num_edges;
                attr_val = [[] for i in range(len(self.wire_edges))];
                for i,ei in enumerate(edge_index_map):
                    attr_val[ei].append(val[
                        i*per_edge_size:(i+1)*per_edge_size]);
                attr_val = [np.mean(val, axis=0) for val in attr_val];
                self.wire_attributes[key] = np.array(attr_val);


                #attr_val = np.zeros(len(self.wire_edges));
                #for i,ei in enumerate(edge_index_map):
                #    attr_val[ei] = val[i];
                #self.wire_attributes[key] = attr_val;

    @timethis
    def __center_at_origin(self):
        bbox_min = np.amin(self.wire_vertices, axis=0);
        bbox_max = np.amax(self.wire_vertices, axis=0);
        bbox_center = 0.5 * (bbox_min + bbox_max);
        self.wire_vertices = self.wire_vertices - bbox_center;

    @timethis
    def __apply_vertex_offset(self):
        if "vertex_offset" in self.wire_attributes:
            offsets = self.wire_attributes["vertex_offset"];
            self.wire_vertices += offsets;

    @timethis
    def __copy_attributes(self, wire_network):
        """ Copy any vertex attributes from single cell to tiled wire network.
        """
        for key,val in self.wire_attributes.iteritems():
            wire_network.attributes.add(key, val);

    @property
    @timethis
    def wire_network(self):
        network = WireNetwork();
        network.load(self.wire_vertices, self.wire_edges);
        self.__copy_attributes(network);
        return network;


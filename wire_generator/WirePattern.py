import numpy as np
import os
import sys

# Update path to import PyMesh
py_mesh_path = os.environ.get("PYMESH_PATH");
if py_mesh_path == None:
    raise ImportError("Please set PYMESH_PATH to the correct lib path.");
sys.path.append(os.path.join(py_mesh_path, "lib"));
sys.path.append(os.path.join(py_mesh_path, "swig"));
import PyMesh

from WireNetwork import WireNetwork
from timethis import timethis

class WirePattern(object):
    def __init__(self, vertices, edges):
        self.pattern_vertices = np.array(vertices);
        self.pattern_edges = np.array(edges, dtype=int);
        self.pattern_bbox_min = np.amin(vertices, axis=0);
        self.pattern_bbox_max = np.amax(vertices, axis=0);
        self.pattern_bbox_size = self.pattern_bbox_max - self.pattern_bbox_min;

    @timethis
    def tile(self, reps):
        self.wire_vertices = np.zeros((0, 3));
        self.wire_edges = np.zeros((0, 2), dtype=int);
        for i in range(reps[0]):
            x_inc = self.pattern_bbox_size[0] * i;
            for j in range(reps[1]):
                y_inc = self.pattern_bbox_size[1] * j;
                for k in range(reps[2]):
                    z_inc = self.pattern_bbox_size[2] * k;

                    base_idx = self.wire_vertices.shape[0];
                    vertices = self.pattern_vertices +\
                            np.array([x_inc, y_inc, z_inc]);
                    self.wire_vertices = np.vstack(
                            (self.wire_vertices, vertices));

                    edges = self.pattern_edges + base_idx;
                    self.wire_edges = np.vstack(
                            (self.wire_edges, edges));
        self.__remove_duplicated_vertices();
        self.__remove_duplicated_edges();

    @timethis
    def __remove_duplicated_vertices(self):
        cell_size = np.amin(self.pattern_bbox_size) * 0.1;
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

    @timethis
    def __remove_duplicated_edges(self):
        edges = set([tuple(sorted(edge)) for edge in self.wire_edges]);
        self.wire_edges = np.array(list(edges), dtype=int);

    @property
    @timethis
    def wire_network(self):
        network = WireNetwork();
        network.load(self.wire_vertices, self.wire_edges);
        return network;


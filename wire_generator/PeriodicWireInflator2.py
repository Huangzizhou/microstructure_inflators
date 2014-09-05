import numpy as np

from WireInflator import WireInflator
from WirePattern import WirePattern

import PyCSG

class PeriodicWireInflator(WireInflator):
    def __init__(self, wire_network):
        super(PeriodicWireInflator, self).__init__(wire_network);
        self.__generate_phantom_wire_network();

        # use phantom wire as input for inflator
        self.original_wire_network = self.wire_network;
        self.wire_network = self.phantom_wire_network;

    def __generate_phantom_wire_network(self):
        wire_pattern = WirePattern();
        wire_pattern.set_single_cell_from_wire_network(self.wire_network);
        wire_pattern.tile([3, 3, 3]);
        phantom_wire_network = wire_pattern.wire_network;

        phantom_wire_network.translate(
                self.wire_network.bbox_center - phantom_wire_network.bbox_center);
        self.phantom_wire_network = phantom_wire_network;

    def inflate(self, clean_up=True):
        if not clean_up:
            raise NotImplementedError("Clean up is required for periodic wires");
        super(PeriodicWireInflator, self).inflate(clean_up);
        self._clip_with_bbox();
        self._enforce_periodic_connectivity();

    def _clip_with_bbox(self):
        assert(self.original_wire_network.dim == 3);
        bbox_min, bbox_max = self.original_wire_network.bbox;
        bbox_vertices, bbox_faces = self._generate_box_mesh(bbox_min, bbox_max);

        csg_engine = PyCSG.CSGEngine.create("cork");
        csg_engine.set_mesh_1(bbox_vertices, bbox_faces);
        csg_engine.set_mesh_2(self.mesh_vertices, self.mesh_faces);
        csg_engine.compute_intersection();

        self.mesh_vertices = csg_engine.get_vertices();
        self.mesh_faces    = csg_engine.get_faces();

    def _generate_box_mesh(self, bbox_min, bbox_max):
        vertices = np.array([
            [bbox_min[0], bbox_min[1], bbox_min[2]],
            [bbox_max[0], bbox_min[1], bbox_min[2]],
            [bbox_max[0], bbox_max[1], bbox_min[2]],
            [bbox_min[0], bbox_max[1], bbox_min[2]],
            [bbox_min[0], bbox_min[1], bbox_max[2]],
            [bbox_max[0], bbox_min[1], bbox_max[2]],
            [bbox_max[0], bbox_max[1], bbox_max[2]],
            [bbox_min[0], bbox_max[1], bbox_max[2]] ], order="C");
        faces = np.array([
            [0, 3, 1],
            [1, 3, 2],
            [4, 5, 7],
            [7, 5, 6],
            [5, 1, 2],
            [5, 2, 6],
            [4, 7, 3],
            [4, 3, 0],
            [0, 1, 5],
            [0, 5, 4],
            [2, 3, 6],
            [3, 7, 6] ], dtype=int, order="C");
        return vertices, faces;

    def _enforce_periodic_connectivity(self):
        self._enforce_single_axis_periodicity(0);
        self._clean_up();
        self._enforce_single_axis_periodicity(1);
        self._clean_up();
        self._enforce_single_axis_periodicity(2);
        self._clean_up();

    def _enforce_single_axis_periodicity(self, axis):
        tol = 1e-3;
        bbox_min, bbox_max = self.original_wire_network.bbox;

        vertices = self.mesh_vertices;
        faces = self.mesh_faces;
        num_vertices = len(self.mesh_vertices);

        v_on_min_boundary = vertices[:, axis] <= bbox_min[axis] + tol;
        v_on_max_boundary = vertices[:, axis] >= bbox_max[axis] - tol;

        offset = np.zeros(self.wire_network.dim);
        offset[axis] = bbox_max[axis] - bbox_min[axis];

        f_on_min_boundary = np.all(v_on_min_boundary[faces], axis=1);
        f_on_max_boundary = np.all(v_on_max_boundary[faces], axis=1);

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

    def write_debug_wires(self, filename, vertices, edges):
        writer = WireWriter(filename);
        writer.write(vertices, edges);

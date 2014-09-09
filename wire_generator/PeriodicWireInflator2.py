import numpy as np
from numpy.linalg import norm
from math import ceil,log

from WireInflator import WireInflator
from WirePattern import WirePattern
from Triangulation import Triangulation

import PyCSG
import PyMeshUtils
import PyMesh

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
        self._subdivide(2);

    def _clip_with_bbox(self):
        assert(self.original_wire_network.dim == 3);
        bbox_min, bbox_max = self.original_wire_network.bbox;
        bbox_vertices, bbox_faces = self._generate_box_mesh(bbox_min, bbox_max);
        #bbox_vertices, bbox_faces = self._subdivide_bbox(bbox_vertices, bbox_faces);

        csg_engine = PyCSG.CSGEngine.create("cork");
        csg_engine.set_mesh_1(bbox_vertices, bbox_faces);
        csg_engine.set_mesh_2(self.mesh_vertices, self.mesh_faces);
        csg_engine.compute_intersection();

        self.mesh_vertices = csg_engine.get_vertices();
        self.mesh_faces    = csg_engine.get_faces();

        self._clean_up();
        self._post_clip_processing(0);
        self._post_clip_processing(1);
        self._post_clip_processing(2);

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

    def _subdivide_bbox(self, vertices, faces):
        bbox_min = np.amin(vertices, axis=0);
        bbox_max = np.amax(vertices, axis=0);
        order = int(ceil(log(
            np.amax(bbox_max-bbox_min) / np.amin(self.thickness), 2)));

        subdiv = PyMeshUtils.Subdivision.create("simple");
        subdiv.subdivide(vertices, faces, order+1);
        return subdiv.get_vertices(), subdiv.get_faces();

    def _enforce_periodic_connectivity(self):
        self._enforce_single_axis_periodicity(0);
        self._clean_up();
        self._enforce_single_axis_periodicity(1);
        self._clean_up();
        self._enforce_single_axis_periodicity(2);
        self._clean_up();

    def _enforce_single_axis_periodicity(self, axis):
        tol = 1e-3;
        axis_dir = np.zeros(3);
        axis_dir[axis] = 1;
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

        added_vertices, added_faces = self._retriangulate(
                vertices, faces[f_on_min_boundary], axis_dir);

        f_on_boundary = np.logical_or(f_on_min_boundary, f_on_max_boundary);
        faces = faces[np.logical_not(f_on_boundary)];

        vertices = np.vstack([vertices, added_vertices, added_vertices + offset]);
        faces = np.vstack([faces, 
            added_faces + num_vertices,
            added_faces[:,[1, 0, 2]] + num_vertices + len(added_vertices) ]);

        #added_vertices = vertices[v_on_min_boundary] + offset;
        #from_index = np.arange(num_vertices, dtype=int)[v_on_min_boundary];
        #to_index = np.arange(len(added_vertices), dtype=int) + num_vertices;
        #vertex_map = {i:j for i,j in zip(from_index, to_index)};
        #index_lookup = lambda i: vertex_map[i];
        #added_faces = [map(index_lookup, face) for face in faces[f_on_min_boundary]];
        #added_faces = np.array(added_faces)[:,[1,0,2]];

        #vertices = np.vstack((vertices, added_vertices));
        #faces = faces[np.logical_not(f_on_max_boundary)];
        #faces = np.vstack((faces, added_faces));

        self.mesh_vertices = vertices;
        self.mesh_faces = faces;

    def _retriangulate(self, vertices, faces, proj_dir):
        """ Retriangulate the polygon formed by the faces without introducing
        new vertices on the boundary.  In particular, all faces should be
        coplanar.
        """
        tri = Triangulation(vertices, faces, proj_dir);
        tri.triangulate(0.1);
        return tri.vertices, tri.faces;

    def _post_clip_processing(self, axis):
        tol = 1e-3;
        bbox_min, bbox_max = self.original_wire_network.bbox;

        offset = np.zeros(self.original_wire_network.dim);
        offset[axis] = (bbox_max - bbox_min)[axis];

        v_on_min_boundary = self.mesh_vertices[:, axis] <= bbox_min[axis] + tol;
        v_on_max_boundary = self.mesh_vertices[:, axis] >= bbox_max[axis] - tol;

        f_on_min_boundary = np.all(v_on_min_boundary[self.mesh_faces], axis=1);
        f_on_max_boundary = np.all(v_on_max_boundary[self.mesh_faces], axis=1);

        min_bd_vertices, min_bd_edges = self._extract_face_boundary(
                self.mesh_vertices, self.mesh_faces[f_on_min_boundary]);
        max_bd_vertices, max_bd_edges = self._extract_face_boundary(
                self.mesh_vertices, self.mesh_faces[f_on_max_boundary]);

        boundary_marker = np.zeros(len(self.mesh_vertices), dtype=int);
        boundary_marker[min_bd_vertices] = -1;
        boundary_marker[max_bd_vertices] = 1;
        self.write_debug_mesh("bd_vertices.msh", boundary_marker);

        min_to_max_map = self._map(
                self.mesh_vertices[min_bd_vertices] + offset,
                self.mesh_vertices[max_bd_vertices]);

        matched = np.zeros(len(self.mesh_vertices), dtype=bool);
        for v,neighbors in zip(min_bd_vertices, min_to_max_map):
            if len(neighbors) == 0: continue;
            elif len(neighbors) == 1:
                other = max_bd_vertices[neighbors[0]];
                matched[v] = True;
                matched[other] = True;
            else:
                source = self.mesh_vertices[v];
                target = self.mesh_vertices[
                        max_bd_vertices[neighbors] ];
                dist = norm(target - source, axis=1);
                idx = np.argmin(dist);
                other = max_bd_vertices[neighbors[idx]];
                matched[v] = True;
                matched[other] = True;

        self.write_debug_mesh("matched_bd.msh", matched.astype(float));

        while not np.all(matched[min_bd_vertices]):
            for edge in min_bd_edges:
                if matched[edge[0]] and matched[edge[1]]: continue;
                if not matched[edge[0]] and matched[edge[1]]:
                    self.mesh_vertices[edge[0]] = self.mesh_vertices[edge[1]];
                    matched[edge[0]] = True;
                elif matched[edge[0]] and not matched[edge[1]]:
                    self.mesh_vertices[edge[1]] = self.mesh_vertices[edge[0]];
                    matched[edge[1]] = True;

        while not np.all(matched[max_bd_vertices]):
            for edge in max_bd_edges:
                if matched[edge[0]] and matched[edge[1]]: continue;
                if not matched[edge[0]] and matched[edge[1]]:
                    self.mesh_vertices[edge[0]] = self.mesh_vertices[edge[1]];
                    matched[edge[0]] = True;
                elif matched[edge[0]] and not matched[edge[1]]:
                    self.mesh_vertices[edge[1]] = self.mesh_vertices[edge[0]];
                    matched[edge[1]] = True;

        self._clean_up();

    def _extract_face_boundary(self, vertices, faces):
        mesh = self._form_mesh(vertices, faces);
        bd_extractor = PyMeshUtils.Boundary.extract_surface_boundary(mesh);
        return  bd_extractor.get_boundary_nodes().ravel(),\
                bd_extractor.get_boundaries();

    def _form_mesh(self, vertices, faces):
        voxels = np.array([]);
        factory = PyMesh.MeshFactory();
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                voxels,
                self.original_wire_network.dim, 3, 0);
        return factory.create();

    def _map(self, from_pts, to_pts):
        dim = from_pts.shape[1];
        num_pts = len(to_pts);
        grid = PyMesh.HashGrid.create(1e-3, dim);
        grid.insert_multiple(np.arange(num_pts, dtype=int), to_pts);

        nearby_pts = [ grid.get_items_near_point(p).ravel() for p in from_pts ];
        return nearby_pts;

    def write_debug_wires(self, filename, vertices, edges):
        writer = WireWriter(filename);
        writer.write(vertices, edges);

    def write_debug_mesh(self, filename, field=None):
        mesh = self._form_mesh(self.mesh_vertices, self.mesh_faces);
        if field is not None:
            mesh.add_attribute("debug");
            mesh.set_attribute("debug", field);

        writer = PyMesh.MeshWriter.create_writer(filename);
        writer.with_attribute("debug");
        writer.write_mesh(mesh);


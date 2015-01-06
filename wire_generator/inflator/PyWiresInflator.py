import numpy as np

import PyWires
import PyMesh

class PyWiresInflator(object):
    def __init__(self, wire_network, parameters, periodic=False):
        self.wire_network = wire_network;
        self.periodic = periodic;
        self.parameters = parameters;

    def inflate(self, clean_up=True, subdivide_order=1,
            subdivide_method="simple",
            rel_geometry_correction=None,
            abs_geometry_correction=None,
            geometry_correction_cap=None):
        wires = self.wire_network.raw_wires;

        if self.periodic:
            inflator = PyWires.InflatorEngine.create_parametric(wires,
                    self.parameters.raw_parameters);
        else:
            inflator = PyWires.InflatorEngine.create("simple", wires);

        thickness = wires.get_attribute("thickness").ravel();
        inflator.set_thickness_type(
                self.parameters.raw_parameters.get_thickness_type());
        inflator.set_thickness(thickness);

        inflator.with_refinement(subdivide_method, subdivide_order);
        if (rel_geometry_correction is not None):
            inflator.with_rel_geometry_correction(np.array(rel_geometry_correction));
        if (abs_geometry_correction is not None):
            inflator.with_abs_geometry_correction(np.array(abs_geometry_correction));
        if (geometry_correction_cap is not None):
            inflator.set_geometry_correction_cap(np.array(geometry_correction_cap));
        inflator.inflate();
        self.mesh_vertices = inflator.get_vertices();
        self.mesh_faces = inflator.get_faces();
        self.source_wire_id = inflator.get_face_sources();

        if self.periodic:
            self.wire_network.compute_symmetry_orbits();

    @property
    def mesh(self):
        factory = PyMesh.MeshFactory();
        factory.load_data(
                self.mesh_vertices.ravel(order="C"),
                self.mesh_faces.ravel(order="C"),
                np.zeros(0),
                self.wire_network.dim, 3, 4);
        mesh = factory.create_shared();
        mesh.add_attribute("source_wire_id");
        mesh.set_attribute("source_wire_id", self.source_wire_id);

        if self.wire_network.has_attribute("vertex_symmetry_orbit"):
            indices = self.wire_network.get_attribute("vertex_symmetry_orbit").ravel();
            source_id_mask = self.source_wire_id > 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("orthotropic_vertex_orbit");
            mesh.set_attribute("orthotropic_vertex_orbit", source_index);

        if self.wire_network.has_attribute("vertex_cubic_symmetry_orbit"):
            indices = self.wire_network.get_attribute("vertex_cubic_symmetry_orbit").ravel();
            source_id_mask = self.source_wire_id > 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("isotropic_vertex_orbit");
            mesh.set_attribute("isotropic_vertex_orbit", source_index);

        if self.wire_network.has_attribute("edge_symmetry_orbit"):
            indices = self.wire_network.get_attribute("edge_symmetry_orbit").ravel();
            source_id_mask = self.source_wire_id < 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[-self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("orthotropic_edge_orbit");
            mesh.set_attribute("orthotropic_edge_orbit", source_index);

        if self.wire_network.has_attribute("edge_cubic_symmetry_orbit"):
            indices = self.wire_network.get_attribute("edge_cubic_symmetry_orbit").ravel();
            source_id_mask = self.source_wire_id < 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[-self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("isotropic_edge_orbit");
            mesh.set_attribute("isotropic_edge_orbit", source_index);

        return mesh;


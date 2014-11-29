import numpy as np

import PyWires
import PyMesh

class PyWiresInflator(object):
    def __init__(self, wire_network, parameters, periodic=False):
        self.wire_network = wire_network;
        self.parameters = parameters;
        self.periodic = periodic;

    def inflate(self, clean_up=True, subdivide_order=1, subdivide_method="simple"):
        wires = PyWires.WireNetwork.create_raw(
                self.wire_network.vertices,
                self.wire_network.edges);
        if self.periodic:
            manager = PyWires.ParameterManager.create(wires);
            inflator = PyWires.InflatorEngine.create_parametric(wires, manager);
        else:
            inflator = PyWires.InflatorEngine.create("simple", wires);

        if "vertex_thickness" in self.wire_network.attributes:
            thickness = self.wire_network.attributes["vertex_thickness"];
            inflator.set_thickness_type(PyWires.InflatorEngine.PER_VERTEX);
            inflator.set_thickness(thickness);
        elif "edge_thickness" in self.wire_network.attributes:
            thickness = self.wire_network.attributes["edge_thickness"];
            inflator.set_thickness_type(PyWires.InflatorEngine.PER_EDGE);
            inflator.set_thickness(thickness);

        inflator.with_refinement(subdivide_method, subdivide_order);
        inflator.inflate();
        self.mesh_vertices = inflator.get_vertices();
        self.mesh_faces = inflator.get_faces();
        self.source_wire_id = inflator.get_face_sources();

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

        if "symmetry_vertex_orbit" in self.wire_network.attributes:
            indices = self.wire_network.attributes["symmetry_vertex_orbit"].ravel();
            source_id_mask = self.source_wire_id > 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("vertex_orbit");
            mesh.set_attribute("vertex_orbit", source_index);

        if "symmetry_edge_orbit" in self.wire_network.attributes:
            indices = self.wire_network.attributes["symmetry_edge_orbit"].ravel();
            source_id_mask = self.source_wire_id < 0;
            source_index = np.zeros_like(self.source_wire_id);
            source_index[np.logical_not(source_id_mask)] = -1;
            source_index[source_id_mask] =\
                    indices[-self.source_wire_id[source_id_mask] - 1];
            mesh.add_attribute("edge_orbit");
            mesh.set_attribute("edge_orbit", source_index);

        return mesh;


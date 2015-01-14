from core.WireNetwork import WireNetwork
import numpy as np

import PyWires
import PyMesh

from InflatorFacade import InflatorFacade

class InflatorFacadeMixedPattern(InflatorFacade):
    def __init__(self, wire_networks):
        assert(len(wire_networks) > 0);
        self.wire_networks = [wires.raw_wires for wires in wire_networks];

    def inflate_with_mixed_patterns(self, mesh, options):
        assert(mesh.has_attribute("pattern_id"));
        tiled_wire_network = self.__tile(mesh, options);
        self.__apply_vertex_offset(tiled_wire_network);
        mesh = self.__inflate(tiled_wire_network, options);
        return mesh;

    def __tile(self, mesh, options):
        tiler = PyWires.WireTiler(self.wire_networks[0]);
        tiled_raw_wires = tiler.tile_with_mixed_patterns(
                self.wire_networks, mesh,
                options.get("thickness_type", "vertex") == "vertex",
                options.get("dof_type", "isotropic") == "isotropic");

        tiled_wire_network = WireNetwork();
        tiled_wire_network.load_from_raw(tiled_raw_wires);
        return tiled_wire_network;

    def __inflate(self, wire_network, options):
        if options.get("trim", False):
            wire_network.trim();

        if options.get("thickness_type", "vertex") == "vertex":
            thickness_type = PyWires.VERTEX;
        else:
            thickness_type = PyWires.EDGE;

        inflator = PyWires.InflatorEngine.create("simple",
                wire_network.raw_wires);
        inflator.set_thickness_type(thickness_type);
        inflator.set_thickness(wire_network.get_attribute("thickness").ravel());
        inflator.with_refinement(options["subdiv_method"], options["subdiv"]);
        self.__apply_geometry_correction(inflator, options);

        inflator.inflate();

        self.mesh_vertices = inflator.get_vertices();
        self.mesh_faces = inflator.get_faces();
        return self.mesh;

    def __apply_geometry_correction(self, inflator, options):
        rel_geometry_correction = options.get("rel_geometry_correction");
        abs_geometry_correction = options.get("abs_geometry_correction");
        geometry_correction_cap = options.get("geometry_correction_cap");
        geometry_spread = options.get("geometry_spread");
        geometry_correction_lookup = options.get("geometry_correction_lookup");

        if (rel_geometry_correction is not None):
            inflator.with_rel_geometry_correction(np.array(rel_geometry_correction));
        if (abs_geometry_correction is not None):
            inflator.with_abs_geometry_correction(np.array(abs_geometry_correction));
        if (geometry_correction_cap is not None):
            inflator.set_geometry_correction_cap(geometry_correction_cap);
        if (geometry_spread is not None):
            inflator.set_geometry_spread_constant(geometry_spread);
        if (geometry_correction_lookup is not None):
            inflator.with_geometry_correction_lookup(geometry_correction_lookup);


    def __apply_vertex_offset(self, wire_network):
        vertices = wire_network.vertices;
        offset = wire_network.get_attribute("vertex_offset");
        wire_network.vertices = vertices + offset;

    @property
    def mesh(self):
        dim = self.mesh_vertices.shape[1];
        vertex_per_face = self.mesh_faces.shape[1];
        vertex_per_voxel = 4;
        factory = PyMesh.MeshFactory();
        factory.load_data(
                self.mesh_vertices.ravel(order="C"),
                self.mesh_faces.ravel(order="C"),
                np.zeros(0),
                dim, vertex_per_face, vertex_per_voxel);
        mesh = factory.create_shared();
        return mesh;


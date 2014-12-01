import numpy as np
from PatternParameter import PatternParameter

class VertexOffsetParameter(PatternParameter):
    def __init__(self, wire_network, orbit_id):
        super(VertexOffsetParameter, self).__init__(wire_network);
        self.default_offset = np.zeros(wire_network.dim);
        self.orbit_id = orbit_id;

        #self.__compute_mask();

    def __compute_mask(self):
        raise DeprecationWarning("This method is deprecated.");
        if "symmetry_vertex_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();
        orbits = self.wire_network.attributes["symmetry_vertex_orbit"];
        orbit_vertices = self.wire_network.vertices[orbits == self.orbit_id];

        bbox_min, bbox_max = self.wire_network.bbox;
        orbit_bbox_min = np.amin(orbit_vertices, axis=0);
        orbit_bbox_max = np.amax(orbit_vertices, axis=0);

        self.__dof_mask = np.logical_or(np.logical_or(
                orbit_bbox_min <= bbox_min, orbit_bbox_max >= bbox_max),
                orbit_bbox_max == orbit_bbox_min);

    def evaluate(self, **kwargs):
        if hasattr(self, "formula"):
            offset = [];
            for component in self.formula:
                if isinstance(component, (str, unicode)):
                    offset.append(eval(component.format(**kwargs)));
                elif isinstance(component, float):
                    offset.append(component);
                else:
                    raise RuntimeError(
                            "Unable to evaluate formula {}".format(
                        self.formula));
            return np.array(offset);
        else:
            return self.default_offset;

    @property
    def dof_mask(self):
        return self.__dof_mask;

    @property
    def names(self):
        dim = self.wire_network.dim;
        coordinate_names = ["x", "y", "z"];
        return ["vertex_orbit_{}_offset_{}".format(
            self.orbit_id, coordinate_names[i]) for i in range(dim)];

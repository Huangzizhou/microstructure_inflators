class WireAttributes(object):
    def __init__(self, wire_network):
        self.__attributes = {};
        self.__wire_network = wire_network;

    def add(self, name, value=None):
        if value is None:
            value = self.__compute_value_from_wire_network(name);
        self.__attributes[name] = value;

    def __compute_value_from_wire_network(self, name):
        if name == "symmetry_vertex_orbit":
            from WireSymmetryVertexOrbitAttribute import WireSymmetryVertexOrbitAttribute
            attr = WireSymmetryVertexOrbitAttribute();
        elif name == "symmetry_edge_orbit":
            from WireSymmetryEdgeOrbitAttribute import WireSymmetryEdgeOrbitAttribute
            attr = WireSymmetryEdgeOrbitAttribute();
        elif name == "isotropic_symmetry_vertex_orbit":
            from WireSymmetryVertexOrbitAttribute import WireSymmetryVertexOrbitAttribute
            attr = WireSymmetryVertexOrbitAttribute();
            attr.set_symmetry_type("isotropic");
        else:
            raise NotImplementedError(
                    "Unknow wire network attribute: {}".format(name));
        attr.compute(self.__wire_network);
        value = attr.value;
        return value;

    def __getitem__(self, key):
        return self.__attributes[key];

    def __setitem__(self, key, val):
        self.__attributes[key] = val;

    def __iter__(self):
        return self.__attributes.__iter__();

    def __delitem__(self, key):
        del self.__attributes[key];

    def __contains__(self, key):
        return key in self.__attributes;

    def keys(self):
        return self.__attributes.keys();

    def iteritems(self):
        return self.__attributes.iteritems();

    def __len__(self):
        return len(self.__attributes);

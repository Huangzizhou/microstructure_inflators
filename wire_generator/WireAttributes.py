class WireAttributes(object):
    def __init__(self, wire_network):
        self.__attributes = {};
        self.__wire_network = wire_network;

    def add(self, name, value=None):
        if name == "symmetry_orbit":
            from WireSymmetryOrbitAttribute import WireSymmetryOrbitAttribute
            attr = WireSymmetryOrbitAttribute();
        else:
            raise NotImplementedError(
                    "Unknow wire network attribute: {}".format(name));

        if value is None:
            attr.compute(self.__wire_network);
        else:
            attr.value = value;
        self.attributes[name] = attr;

    @property
    def attributes(self):
        return self.__attributes;

    def __getitem__(self, key):
        return self.__attributes[key];

    def __setitem__(self, key, val):
        self.__attributes[key] = val;

    def __delitem__(self, key):
        del self.__attributes[key];

    def __contains__(self, key):
        return key in self.__attributes;

    def keys(self):
        return self.__attributes.keys();


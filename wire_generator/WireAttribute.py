
class WireAttribute(object):
    def compute(self, wire_network):
        raise NotImplementedError("This method is abstract");

    @property
    def value(self):
        raise NotImplementedError("This method is abstract");

    @value.setter
    def value(self):
        raise NotImplementedError("This method is abstract");

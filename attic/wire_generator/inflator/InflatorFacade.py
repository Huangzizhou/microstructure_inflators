from WireTiler import WireTiler

class InflatorFacade(object):
    @classmethod
    def create(cls, wire_network, parameters):
        dim = wire_network.dim;
        if dim == 2:
            import InflatorFacade2D;
            return InflatorFacade2D.InflatorFacade2D(wire_network, parameters);
        elif dim == 3:
            import InflatorFacade3D;
            return InflatorFacade3D.InflatorFacade3D(wire_network, parameters);
        else:
            raise NotImplementedError("Unsupported dimention: {}".format(dim));

    @classmethod
    def create_mixed(cls, wire_networks):
        import InflatorFacadeMixedPattern;
        return InflatorFacadeMixedPattern.InflatorFacadeMixedPattern(
                wire_networks);

    def __init__(self, wire_network, parameters):
        self.unit_pattern = wire_network;
        self.parameters = parameters;

    def inflate_with_guide_box(self, bbox_min, bbox_max, repetitions, options):
        """ Valid options:
        {
            "trim": bool,
            "periodic": bool
        }
        """
        raise NotImplementedError("This method is abstract");

    def inflate_with_guide_mesh(self, mesh, options):
        """ Valid options:
        {
            "trim": bool
        }
        """
        raise NotImplementedError("This method is abstract");

    def inflate_with_mixed_patterns(self, mesh, options):
        """ Valid options:
        {
            "trim": bool
        }
        """
        raise NotImplementedError("This method is abstract");


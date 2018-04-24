from ThicknessParameter import ThicknessParameter

class VertexThicknessParameter(ThicknessParameter):
    @property
    def names(self):
        return ["vertex_orbit_{}_thickness".format(self.orbit_id)];


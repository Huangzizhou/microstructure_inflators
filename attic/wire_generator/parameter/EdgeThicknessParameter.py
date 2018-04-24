from ThicknessParameter import ThicknessParameter

class EdgeThicknessParameter(ThicknessParameter):
    @property
    def names(self):
        return ["edge_orbit_{}_thickness".format(self.orbit_id)];


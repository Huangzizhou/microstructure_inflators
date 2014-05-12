from ggplot import ggplot
from AbstractPlot import AbstractPlot

class Histogram(AbstractPlot):
    def __init__(self, plot, config):
        self.plot = plot;
        self.config = config;

    def process_plot(self):
        """ syntax:
        "x": "column_1",
        "num_bins": 20,
        "fill": "column_2",
        "position": dodge
        """
        column = self.config["x"];
        num_bins = self.config["num_bins"];
        fill = self.config.get("fill", None);
        position = self.config.get("position", None);
        if position not in [None, "stack", "dodge"]:
            raise NotImplementedError("Unknown position specification: {}"\
                    .format(position));
        self.plot.histogram(column, num_bins, fill, position);

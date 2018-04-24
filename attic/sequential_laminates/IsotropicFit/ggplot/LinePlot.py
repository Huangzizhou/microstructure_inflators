from ggplot import ggplot
from AbstractPlot import AbstractPlot

class LinePlot(AbstractPlot):
    def __init__(self, plot, config):
        self.plot = plot;
        self.config = config;

    def process_plot(self):
        """ syntax:
        "x": "column_1",
        "y": "column_2",
        "w": "column_3",
        "group": "group",
        "with_points": bool
        """
        x = self.config.get("x");
        y = self.config.get("y");
        w = self.config.get("w");
        group = self.config.get("group");

        self.plot.line_plot(x, y, w, group);
        if self.config.get("with_points"):
            self.plot.scatter_plot(x, y, w);

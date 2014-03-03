from ggplot import ggplot
from AbstractPlot import AbstractPlot

class ContourPlot(AbstractPlot):
    def __init__(self, plot, config):
        self.plot = plot;
        self.config = config;

    def process_plot(self):
        """ syntax:
        "x": "column_1",
        "y": "column_2",
        "w": "column_3",
        "with_points": bool,
        "num_bins": #
        """
        x = self.config.get("x");
        y = self.config.get("y");
        w = self.config.get("w");
        num_bins = self.config.get("num_bins", 10);

        self.plot.contour_plot(x, y, w, num_bins);
        if self.config.get("with_points"):
            self.plot.scatter_plot(x, y, w);

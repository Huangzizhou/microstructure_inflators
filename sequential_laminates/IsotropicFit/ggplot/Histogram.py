from ggplot import ggplot
from AbstractPlot import AbstractPlot

class Histogram(AbstractPlot):
    def __init__(self, plot, config):
        self.plot = plot;
        self.config = config;

    def process_plot(self):
        """ syntax:
        "x": "column_1",
        "num_bins": 20
        """
        column = self.config["x"];
        num_bins = self.config["num_bins"];
        self.plot.histogram(column, num_bins);

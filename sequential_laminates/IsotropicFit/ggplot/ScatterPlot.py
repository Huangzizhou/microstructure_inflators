from ggplot import ggplot
from AbstractPlot import AbstractPlot

class ScatterPlot(AbstractPlot):
    def __init__(self, plot, config):
        self.plot = plot;
        self.config = config;

    def process_plot(self):
        """ syntax:
        "x": "column_1",
        "y": "column_2",
        "w": "column_3",
        "with_lines": bool
        """
        x = self.config.get("x");
        y = self.config.get("y");
        w = self.config.get("w");

        self.plot.scatter_plot(x, y, w);
        if self.config.get("with_lines"):
            self.plot.line_plot(x, y, w);

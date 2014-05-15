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
        "with_lines": bool,
        "shape": shape_column,
        "point_size": point_size_column
        """
        x = self.config.get("x");
        y = self.config.get("y");
        w = self.config.get("w");
        shape = self.config.get("shape", None);
        point_size = self.config.get("point_size", None);
        group = self.config.get("group");

        self.plot.scatter_plot(x, y, w, shape, point_size);
        if self.config.get("with_lines"):
            self.plot.line_plot(x, y, w, group);

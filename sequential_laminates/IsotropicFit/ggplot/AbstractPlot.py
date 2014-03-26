import os.path
from ggplot import ggplot
from deprecated import deprecated

class AbstractPlot:
    def __init__(self, plot, config):
        self.plot = plot;
        self.config = config;

    def generate_r_script(self):
        self.preprocess_plot();
        self.process_plot();
        self.postprocess_plot();

    def preprocess_plot(self):
        self.discretize_column();
        self.add_column();
        self.filter_data();
        self.init_plot();

    def process_plot(self):
        raise NotImplementedError("process_plot() is abstract");

    def postprocess_plot(self):
        self.add_facet_grid();
        self.add_points();
        self.add_lines();
        self.add_texts();
        self.set_range();
        self.log_scale();
        self.set_title();
        self.coord_fixed();
        self.save_plot();

    @deprecated
    def discretize_column(self):
        """ syntax:
        discretize_columns : [
            {"name": "column_1"},
            {"name": "column_2"},
            ...
        ]
        """
        columns = self.config.get("discretize_columns", None);
        if columns is not None:
            for column in columns:
                self.plot.discretize_column(column["name"]);

    @deprecated
    def add_column(self):
        """ syntax:
        add_columns : [
            {"name": "column_1", "formula": "column_0 * 2"},
            {"name": "column_2", "formula": "column_0 + column_1"},
            ...
        ]
        """
        columns = self.config.get("add_columns", None);
        if columns is not None:
            for column in columns:
                self.plot.add_column(column["name"], column["formula"]);

    @deprecated
    def filter_data(self):
        """ syntax:
        add_filters : [
            {"name": "column_1", "min_val": -1, "max_val": 1},
            {"name": "column_2", "min_val":  0, "max_val": 5},
            ...
        ]
        """
        filters = self.config.get("add_filters", None);
        if filters is not None:
            for f in filters:
                self.plot.filter_data(f["name"],
                        min_val = f["min_val"],
                        max_val = f["max_val"]);

    def init_plot(self):
        self.plot.init_plot();

    def add_facet_grid(self):
        """ syntax:
        "facet_1" : "column_1",
        "facet_2" : "column_2",
        """
        facet_1 = self.config.get("facet_1", ".");
        facet_2 = self.config.get("facet_2", ".");
        if facet_1 != "." or facet_2 != ".":
            self.plot.facet_grid(facet_1, facet_2);

    def set_range(self):
        """ syntax:
        "x_min" : -1,
        "x_max" :  1,
        "y_min" :  0,
        "y_max" :  5
        """
        self.plot.set_range(
                x_min = self.config.get("x_min", None),
                x_max = self.config.get("x_max", None),
                y_min = self.config.get("y_min", None),
                y_max = self.config.get("y_max", None));

    def log_scale(self):
        """ syntax:
        "log_scale": "x" | "y" | "xy"
        """
        log_axis = self.config.get("log_scale", "");
        self.plot.log_scale(
                x = "x" in log_axis,
                y = "y" in log_axis);

    def set_title(self):
        """ syntax:
        "title" : "x vs y"
        """
        title = self.config.get("title", None);
        if title is not None:
            self.plot.title(title);

    def add_points(self):
        """ syntax:
        "add_points" : [
            {"x": 0, "y": 1, "color":"green", "size":3},
            {"x": 1, "y": 1, "color":"red"},
            {"x": 2, "y": 5}
        ]
        """
        points = self.config.get("add_points", None);
        if points is not None:
            for p in points:
                self.plot.draw_point(**p);

    def add_lines(self):
        self.add_hline();
        self.add_vline();
        self.add_line();

    def add_hline(self):
        """ syntax:
        "hline": y
        or
        "hline": [y1, y2, ...]
        """
        y = self.config.get("hline", None);
        if y is not None:
            if not isinstance(y, list):
                y = [y];
            for y_intercept in y:
                self.plot.draw_hline(y_intercept);

    def add_vline(self):
        """ syntax:
        "vline": x
        or
        "vline": [x1, x2, ...]
        """
        x = self.config.get("vline", None);
        if x is not None:
            if not isinstance(x, list):
                x = [x];
            for x_intercept in x:
                self.plot.draw_vline(x_intercept);

    def add_line(self):
        """ syntax:
        "add_lines": [
            {
                "intercept": intercept,
                "slope": slope
            },
            ...
        ]
        """
        line_settings = self.config.get("add_lines", None);
        if line_settings is not None:
            for setting in line_settings:
                intercept = setting["intercept"];
                slope = setting["slope"];
                self.plot.draw_line(intercept, slope);

    def add_texts(self):
        """ syntax:
        "add_texts" : [
            {"x": 0, "y": 0, "text": "origin", "hjust": 0, "vjust": 1},
            {"x": 1, "y": 1, "text": "(1, 1)"}
        ]
        """
        texts = self.config.get("add_texts", None);
        if texts is not None:
            for t in texts:
                self.plot.draw_text(**t);

    def coord_fixed(self):
        """ syntax:
        "coord_fixed": ratio
        """
        ratio = self.config.get("coord_fixed", None);
        if (ratio is not None):
            self.plot.coord_fixed(ratio);

    def save_plot(self):
        """ syntax:
        "out_name": "plot_1.png",
        "width": 10,
        "height": 6
        """
        args = {}
        out_name = self.config["out_name"];
        if not os.path.isabs(out_name):
            out_name = os.path.join(self.config["config_dir"], out_name);
        args["out_name"] = out_name;
        if "width" in self.config:
            args["width"] = self.config["width"];
        if "height" in self.config:
            args["height"] = self.config["height"];
        self.plot.save_plot(**args);




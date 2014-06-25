from ggplot import ggplot

class DataFrame:
    """ Generator of R methods for data manipulation
    Such manipulation would change the data frame and thus influence all plots.
    """

    def __init__(self, plots, config):
        self.plots = plots;
        self.config = config;

    def prepare(self):
        self.add_column();
        self.discretize_column();
        self.as_numeric();
        self.filter_data();

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
                self.plots.discretize_column(column["name"]);

    def as_numeric(self):
        """ syntax:
        as_numeric: [
            {"name": "column_1"},
            {"name": "column_2"},
            ...
        ]
        """
        columns = self.config.get("as_numeric", None);
        if columns is not None:
            for column in columns:
                self.plots.as_numeric(column["name"]);

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
                self.plots.add_column(column["name"], column["formula"]);

    def filter_data(self):
        """ syntax:
        add_filters : [
            {"name": "column_1", "min_val": -1, "max_val": 1},
            {"name": "column_2", "min_val":  0, "max_val": 5},
            ...
        ]
        or 
        add_filters : [
            {"condition": "col1 < col2 and col1 > 0.0"},
            ...
        ]
        """
        filters = self.config.get("add_filters", None);
        if filters is not None:
            for f in filters:
                if "condition" not in f:
                    self.plots.filter_data(f["name"],
                            min_val = f["min_val"],
                            max_val = f["max_val"]);
                else:
                    self.plots.filter_data_with_condition(
                            f["condition"]);



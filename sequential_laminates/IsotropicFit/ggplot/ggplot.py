from csv import DictReader
from mako.template import Template
from mako.runtime import Context
import os.path
from subprocess import check_call
import sys

class ggplot:
    def __init__(self):
        src_dir = sys.path[0]
        self.template = Template(
                filename=os.path.join(src_dir, "ggplot.mako"));
        self.r_script = "";

    def set_data(self, csv_file):
        self.r_script += self.template.get_def("header").render(
                csv_file=csv_file);

    def init_plot(self):
        self.r_script += self.template.get_def("init_plot").render();

    def filter_data(self, col_name, min_val, max_val):
        self.r_script += self.template.get_def("filter_data").render(
                field_name = col_name,
                lower_bound = min_val,
                upper_bound = max_val);

    def filter_data_with_condition(self, condition):
        self.r_script += self.template.get_def(
                "filter_data_with_condition").render(
                        condition = condition);

    def add_column(self, col_name, formula):
        self.r_script += self.template.get_def("add_column").render(
                new_column = col_name, formula = formula);

    def discretize_column(self, col_name):
        self.r_script += self.template.get_def("discretize_column").render(
                field_name = col_name);

    def scatter_plot(self, x_col, y_col, w_col, shape=None, point_size=None, discrete=False):
        self.r_script += self.template.get_def("scatter_plot").render(
                x_col = x_col, y_col = y_col, w_col = w_col, discrete=discrete,
                shape = shape, point_size = point_size);

    def line_plot(self, x_col, y_col, w_col, group=None):
        self.r_script += self.template.get_def("line_plot").render(
                x_col = x_col, y_col = y_col, w_col = w_col, group=group);

    def contour_plot(self, x_col, y_col, w_col, num_bins):
        self.r_script += self.template.get_def("contour_plot").render(
                x_col = x_col, y_col = y_col, w_col = w_col, num_bins = num_bins);

    def title(self, title_text):
        self.r_script += self.template.get_def("title").render(
                title_text = title_text);

    def histogram(self, w_col, num_bins=20, fill=None, position=None):
        self.r_script += self.template.get_def("histogram").render(
                w_col = w_col, num_bins = num_bins, fill = fill,
                position = position);

    def set_range(self, x_min=None, x_max=None, y_min=None, y_max=None):
        self.r_script += self.template.get_def("set_range").render(
                x_col_min = x_min, x_col_max = x_max,
                y_col_min = y_min, y_col_max = y_max);

    def log_scale(self, x=True, y=True):
        self.r_script += self.template.get_def("log_scale").render(
                x=x, y=y);

    def facet_grid(self, facet_1, facet_2):
        self.r_script += self.template.get_def("facet_grid").render(
                facet_1 = facet_1, facet_2 = facet_2);

    def facet_wrap(self, facet):
        self.r_script += self.template.get_def("facet_wrap").render(
                facet = facet);

    def smooth(self, x_col, y_col, method="lm", group=None):
        self.r_script += self.template.get_def("smooth").render(
                x_col = x_col, y_col = y_col, method=method, group=group);

    def save_plot(self, out_name, width=10, height=6):
        self.r_script += self.template.get_def("save_plot").render(
                out_name = out_name, width = width, height = height);

    def draw_point(self, x, y, color="green", size=3):
        self.r_script += self.template.get_def("draw_point").render(
                x=x, y=y, color=color, size=size);

    def draw_hline(self, y):
        self.r_script += self.template.get_def("draw_hline").render(y=y);

    def draw_vline(self, x):
        self.r_script += self.template.get_def("draw_vline").render(x=x);

    def draw_line(self, intercept, slope):
        self.r_script += self.template.get_def("draw_line").render(
                intercept=intercept, slope = slope);

    def draw_text(self, x, y, text, hjust=0, vjust=1):
        self.r_script += self.template.get_def("draw_text").render(
                x=x, y=y, text=text, hjust=hjust, vjust=vjust);

    def coord_fixed(self, ratio=1.0):
        self.r_script += self.template.get_def("coord_fixed").render(
                ratio=ratio);

    def backup_data(self):
        self.r_script += self.template.get_def("backup_data").render();

    def restore_data(self):
        self.r_script += self.template.get_def("restore_data").render();

    def plot(self, r_file):
        print(r_file);
        with open(r_file, 'w') as fout:
            fout.write(self.r_script);
        command = "Rscript {}".format(r_file);
        check_call(command.split());




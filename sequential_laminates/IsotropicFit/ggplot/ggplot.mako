
<%def name="header()">
library("ggplot2")

raw_data <- read.csv("${csv_file}");
names(raw_data) <- sub(" ", "_", names(raw_data));

</%def>

<%def name="backup_data()">
backup = raw_data;
</%def>

<%def name="restore_data()">
raw_data = backup;
</%def>

<%def name="check_nonempty_data()">
cat(nrow(raw_data), "samples left after filtering\n")
if (nrow(raw_data) == 0) {
    print("No data left after filtering");
    q();
}
</%def>

<%def name="filter_data()">
raw_data <- subset(raw_data, ${field_name} < ${upper_bound} & ${field_name} > ${lower_bound});
${check_nonempty_data()}
</%def>

<%def name="filter_data_with_condition()">
raw_data <- subset(raw_data, ${condition});
${check_nonempty_data()}
</%def>

<%def name="add_column()">
raw_data$${new_column} = with(raw_data, ${formula});
</%def>

<%def name="discretize_column()">
raw_data$${field_name} = as.factor(raw_data$${field_name});
</%def>

<%def name="init_plot()">
p <- ggplot(raw_data);
</%def>

<%def name="scatter_plot()">
<%
aes = ["x={}".format(x_col), "y={}".format(y_col)]
if w_col is not None:
    aes.append("color={}".format(w_col));
if shape is not None:
    aes.append("shape={}".format(shape));
if point_size is not None:
    aes.append("size={}".format(point_size));
aes = ",".join(aes);
%>
p <- p + geom_point(aes(${aes}));
% if w_col is not None:
if (is.factor(raw_data$${w_col})) {
    p <- p + scale_colour_brewer(palette="Set1")
} else {
    p <- p + scale_color_gradientn(colours=c("blue", "green", "orange", "red"));
}
% endif
</%def>

<%def name="line_plot()">
<%
clauses = [];
clauses.append("x={}".format(x_col));
clauses.append("y={}".format(y_col));
if w_col is not None:
    clauses.append("color={}".format(w_col));
if group is not None:
    clauses.append("group={}".format(group));
else:
    clauses.append("group=1");
aes = ",".join(clauses);
%>
p <- p + geom_line(aes(${aes}));
% if w_col is not None:
if (is.factor(raw_data$${w_col})) {
    p <- p + scale_colour_brewer(palette="Set1")
} else {
    p <- p + scale_color_gradientn(colours=c("blue", "green", "orange", "red"));
}
% endif
</%def>

<%def name="contour_plot()">
#p <- p + geom_tile(aes(x=${x_col}, y=${y_col}, fill=${w_col}));
#p <- p + scale_fill_gradientn(colours=c("blue", "green", "orange", "red"));
p <- p + stat_contour(aes(x=${x_col}, y=${y_col}, z=${w_col}, color=..level..), bins=${num_bins});
p <- p + labs(color="${w_col}");
#p <- p + scale_color_gradientn(colours=c("blue", "green", "orange", "red"));
</%def>

<%def name="title()">
p <- p + ggtitle("${title_text}");
</%def>

<%def name="histogram()">
bin_width <- (max(raw_data$${w_col}) - min(raw_data$${w_col}))/${num_bins};
<%
clauses = ["x={}".format(w_col)];
if fill is not None:
    clauses.append("fill={}".format(fill));
aes = ",".join(clauses);

opts = ["binwidth=bin_width"];
if position is not None:
    opts.append("position=\"{}\"".format(position));
options = ",".join(opts);
%>
p <- p + geom_histogram(aes(${aes}), ${options});
</%def>


<%def name="set_range()">
% if x_col_min is not None and x_col_max is not None:
p <- p + xlim(${x_col_min}, ${x_col_max});
% endif
% if y_col_min is not None and y_col_max is not None:
p <- p + ylim(${y_col_min}, ${y_col_max});
% endif
</%def>

<%def name="log_scale()">
% if x:
p <- p + scale_x_log10();
% endif
% if y:
p <- p + scale_y_log10();
% endif
</%def>

<%def name="facet_grid()">
p <- p + facet_grid(${facet_1} ~ ${facet_2});
</%def>

<%def name="facet_wrap()">
p <- p + facet_wrap(~${facet});
</%def>

<%def name="smooth()">
<%
clauses = ["x={}".format(x_col), "y={}".format(y_col)]
if group is not None:
    clauses.append("group={}".format(group));
clauses = ",".join(clauses);
%>
p <- p + stat_smooth(aes(${clauses}), method="${method}");
</%def>

<%def name="save_plot()">
ggsave("${out_name}", width=${width}, height=${height});
</%def>

<%def name="draw_ave_points()">
x <- mean(raw_data$${x_col});
y <- mean(raw_data$${y_col});
p <- p + annotate("text", x=x, y=y, label="${label}", hjust=0, vjust=1);
p <- p + annotate("point", x=x, y=y, color="green", size=3);
</%def>

<%def name="draw_point()">
p <- p + annotate("point", x=${x}, y=${y}, color="${color}", size=${size});
</%def>

<%def name="draw_hline()">
p <- p + geom_hline(aes(yintercept=${y}));
</%def>

<%def name="draw_vline()">
p <- p + geom_vline(aes(xintercept=${x}));
</%def>

<%def name="draw_line()">
p <- p + geom_abline(intercept=${intercept}, slope=${slope});
</%def>

<%def name="draw_text()">
p <- p + annotate("text", x=${x}, y=${y}, label="${text}", hjust=${hjust}, vjust=${vjust});
</%def>

<%def name="coord_fixed()">
p <- p + coord_fixed(ratio=${ratio});
</%def>


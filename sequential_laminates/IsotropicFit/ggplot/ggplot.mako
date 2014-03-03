
<%def name="header()">
library("ggplot2")

raw_data <- read.csv("${csv_file}");
names(raw_data) <- sub(" ", "_", names(raw_data));

</%def>

<%def name="filter_data()">
raw_data <- subset(raw_data, ${field_name} < ${upper_bound} & ${field_name} > ${lower_bound});
cat(nrow(raw_data), "samples left after filtering", "${field_name}", "\n")
if (nrow(raw_data) == 0) {
    print("No data left after filtering");
    q();
}
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
% if w_col is not None:
p <- p + geom_point(aes(x=${x_col}, y=${y_col}, color=${w_col}));
% else:
p <- p + geom_point(aes(x=${x_col}, y=${y_col}));
% endif
% if discrete:
#p <- p + scale_color_discrete();
p <- p + scale_colour_brewer(palette="Paired")
% else:
#p <- p + scale_color_gradient(low="blue", high="red");
p <- p + scale_color_gradientn(colours=c("blue", "green", "orange", "red"));
% endif
</%def>

<%def name="line_plot()">
% if w_col is not None:
p <- p + geom_line(aes(x=${x_col}, y=${y_col}, color=${w_col}));
% else:
p <- p + geom_line(aes(x=${x_col}, y=${y_col}));
% endif
</%def>

<%def name="title()">
p <- p + ggtitle("${title_text}");
</%def>

<%def name="histogram()">
bin_width <- (max(raw_data$${w_col}) - min(raw_data$${w_col}))/${num_bins};
p <- p + geom_histogram(aes(x=${w_col}), binwidth=bin_width);
</%def>


<%def name="set_range()">
% if x_col_min is not None and x_col_max is not None:
p <- p + xlim(${x_col_min}, ${x_col_max});
% endif
% if y_col_min is not None and y_col_max is not None:
p <- p + ylim(${y_col_min}, ${y_col_max});
% endif
</%def>

<%def name="facet_grid()">
p <- p + facet_grid(${facet_1} ~ ${facet_2});
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

<%def name="draw_text()">
p <- p + annotate("text", x=${x}, y=${y}, label="${text}", hjust=${hjust}, vjust=${vjust});
</%def>


#!/usr/bin/env perl -w
use strict;

# for my $algo ('simple') {
for my $subdiv (0,1,2) {
# for my $vol (qw(0.001 0.0005 0.00025 0.000125 0.0000625)) {
for (my $numIntervals = 1; $numIntervals <= 20; ++$numIntervals) {
open(GNUPLOT, "| gnuplot") or die();
print GNUPLOT <<PLOT_END;
    set term pngcairo size 800,600;
    set output 'S${subdiv}_n${numIntervals}.png';
    set title 'Simple Subdivision x$subdiv, $numIntervals Intervals';
    unset key;
    set grid;
    set logscale y;
    set xrange [0:30];
    set yrange [1e-9:1000];
    set xlabel 'Iteration';
    set ylabel 'Objective';
PLOT_END
for (my $interval = 0; $interval < $numIntervals; ++$interval) {
for my $uplow ('up','low') {
    my $plotString = ($interval == 0 && $uplow eq 'up') ? "plot " : "";
    my $terminator = ($interval < $numIntervals - 1  || $uplow eq 'up') ? ", \\" : ";";
    print GNUPLOT "$plotString'< grep JS S${subdiv}/num$numIntervals/out_$uplow.$interval.txt | cut -f2' with lines notitle$terminator\n";
}
}
close(GNUPLOT);
}
}
# }
# }

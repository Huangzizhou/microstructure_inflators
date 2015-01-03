#!/usr/bin/env perl -w
use strict;

# for my $algo ('simple') {
for my $subdiv (0,1,2) {
# for my $vol (qw(0.001 0.0005 0.00025 0.000125 0.0000625)) {
for (my $numIntervals = 1; $numIntervals <= 20; ++$numIntervals) {
for (my $interval = 0; $interval < $numIntervals; ++$interval) {
for my $uplow ('up','low') {
open(GNUPLOT, "| gnuplot") or die();
print GNUPLOT <<PLOT_END;
    set term pngcairo size 800,600;
    set output 'S${subdiv}_n${numIntervals}_i${interval}_$uplow.png';
    set title 'Simple Subdivision x$subdiv, Interval $interval/$numIntervals, starting $uplow';
    unset key;
    set grid;
    set logscale y;
    set xlabel 'Iteration';
    set ylabel 'Objective';
    plot '< grep "JS" S${subdiv}/num$numIntervals/out_$uplow.$interval.txt | cut -f2' with lines;
PLOT_END
close(GNUPLOT);
}
}
}
# }
}
# }

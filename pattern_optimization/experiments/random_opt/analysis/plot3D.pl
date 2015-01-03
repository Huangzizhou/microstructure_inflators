#!/usr/bin/env perl -w
use strict;

for my $subdiv (0,1,2,3) {
for my $dist (qw(0.05 0.10 0.20 0.25 0.33 0.50 1.00 16.0)) {
open(GNUPLOT, '| gnuplot');
print GNUPLOT <<END_HEADER;
    set term png size 800,600;
    set output '${subdiv}_${dist}.png';
    set xlabel 'Iteration';
    set ylabel 'Compliance Tensor Fit Objective';
    set logscale y;
    set xrange [0:45];
    set yrange [1e-9:1e1];
    set format y "%.1e";
    set grid;
    unset key;
    set title 'optimization from max-dist $dist away (meshing failures marked in red)';
END_HEADER
for (my $test = 1; $test <= 20; ++$test) {
    my $plotString = $test == 1 ? "plot " : "";
    my $terminator = $test < 20 ? ", \\" : ";";
    print GNUPLOT "$plotString'< grep JS S$subdiv/dist$dist/result.$test.txt | cut -f2' with lines notitle";
    # Mark endpoints, highlighting the failed ones with big circles squares (no Ceres convergence line)
    my $numIters = `grep JS S$subdiv/dist$dist/result.$test.txt | wc -l`;
    $numIters =~ s/\s*//g;
    if (`tail -n1 S$subdiv/dist$dist/result.$test.txt` !~ /Ceres/) {
        print GNUPLOT ", \\\n'< grep JS S$subdiv/dist$dist/result.$test.txt | cut -f2 | tail -n1' using ($numIters - 1):1 with points lc rgb 'black' ps 1.7 pt 7 notitle";
        print GNUPLOT ", \\\n'< grep JS S$subdiv/dist$dist/result.$test.txt | cut -f2 | tail -n1' using ($numIters - 1):1 with points lc rgb 'red' ps 1.5 pt 7 notitle";
    }
    else {
        print GNUPLOT ", \\\n'< grep JS S$subdiv/dist$dist/result.$test.txt | cut -f2 | tail -n1' using ($numIters - 1):1 with points lc rgb 'gray' ps 0.5 pt 7 notitle";
    }
    print GNUPLOT "$terminator\n";
}
close(GNUPLOT);
}
}

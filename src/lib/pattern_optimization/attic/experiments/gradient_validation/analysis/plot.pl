#!/usr/bin/env perl
use strict;
use warnings;

for my $algo ('simple','loop') {
for my $subdiv (0,1,2,3) {
for my $col (3,4) {
for my $vol (qw(0.001 0.0005 0.00025 0.000125 0.0000625)) {
for my $range ('comp7,-0.2,0.2','comp7,-0.02,0.02') {
    my $xrange = $range;
    $xrange =~ s/.*,([^,]+),(.*)/$1:$2/g;
my $type = ($col == 3) ? "objective" : "gradient";
open(GNUPLOT, "| gnuplot") or die();
print GNUPLOT <<PLOT_END;
    set term png size 800,600;
    set output '$type.${range}_${algo}_${subdiv}_$vol.png';
    set title '$algo Subdivision x$subdiv, max vol $vol $range: $type';
    unset key;
    set grid;
    set xrange [$xrange];
    set xlabel 'p';
    set ylabel '$type';
    plot '< grep "^[0-9]" $algo/S$subdiv/vol$vol/result.$range.txt' using 2:$col with lines;
PLOT_END
close(GNUPLOT);
}
}
}
}
}

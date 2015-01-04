#!/usr/bin/env perl
use strict;
use warnings;

print "flippers = [";

for my $algo ('simple','loop') {
    for my $subdiv (0,1,2,3) {
        my $flipperName = "${algo}_x$subdiv.js";
        open(FLIPPER, "> $flipperName");

        my @ranges = ('comp7,-0.2,0.2','comp7,-0.02,0.02');
        my $views = join(', ', map { "'Objective ($_)'" } @ranges) . ", " .
                    join(', ', map { "'Gradient ($_)'" } @ranges);

        print FLIPPER <<END_HEADER;
title = '$algo subdiv $subdiv';
statistics = [];
views = [$views];
frames = [
END_HEADER

        for my $vol (qw(0.001 0.0005 0.00025 0.000125 0.0000625)) {
            my $images = join(', ', map { "'objective.${_}_${algo}_${subdiv}_$vol.png'" } @ranges) . ", " .
                         join(', ', map { "'gradient.${_}_${algo}_${subdiv}_$vol.png'" } @ranges);
            print FLIPPER "{'name': 'Vol $vol', 'image': [$images]},\n";
        }
        print FLIPPER "];\n";
        close(FLIPPER);

        print "    ['$algo x$subdiv', '$flipperName'],\n";
    }
}

print "];\n";

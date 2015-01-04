#!/usr/bin/env perl
use strict;
use warnings;

print "flippers = [";
for my $subdiv (0,1,2,3) {
    open(FLIPPER, "> S$subdiv.js");
    print FLIPPER <<END_HEADER;
title = 'Subdiv $subdiv';
statistics = [];
views = ['Convergence'];
frames = [
END_HEADER
    for my $dist (qw(0.05 0.10 0.20 0.25 0.33 0.50 1.00 16.0)) {
        print FLIPPER "{'name': 'Dist $dist', 'image': ['${subdiv}_$dist.png']},\n"
    }
    print FLIPPER "];\n";
    close FLIPPER;

    print "    ['Sub $subdiv', 'S$subdiv.js'],\n";
}
print "];\n";

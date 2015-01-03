#!/usr/bin/env perl -w
use strict;
print "flippers = [";
for (my $numIntervals = 1; $numIntervals <= 20; ++$numIntervals) {
    for my $subdiv (0,1,2) {
        my $flipperPath = "n$numIntervals.i$subdiv.js";
        open(FLIPPER, ">", $flipperPath);
        print FLIPPER <<END_HEADER;
title = '$numIntervals Intervals, Inflator Subdiv $subdiv';
statistics = [];
views = ['Start at Upper Bound', 'Start at Lower Bound'];
frames = [
END_HEADER
        for (my $interval = 0; $interval < $numIntervals; ++$interval) {
            print FLIPPER "{'name': '$interval', 'image': ['S${subdiv}_n${numIntervals}_i${interval}_up.png', 'S${subdiv}_n${numIntervals}_i${interval}_low.png']},\n";
        }
        print FLIPPER "];";
        close(FLIPPER);
        print "    ['$numIntervals Intervals:Sub $subdiv', '$flipperPath'],\n";
    }
}
print "];\n";

#!/usr/bin/env perl -w
use strict;
print "flippers = [";
for (my $numIntervals = 1; $numIntervals <= 20; ++$numIntervals) {
    for (my $interval = 0; $interval < $numIntervals; ++$interval) {
        my $flipperPath = "n$numIntervals.i$interval.js";
        open(FLIPPER, ">", $flipperPath);
        print FLIPPER <<END_HEADER;
title = 'Interval $interval/$numIntervals';
statistics = [];
views = ['Start at Upper Bound', 'Start at Lower Bound'];
frames = [
END_HEADER
        for my $subdiv (0,1,2) {
            print FLIPPER "{'name': 'S$subdiv', 'image': ['S${subdiv}_n${numIntervals}_i${interval}_up.png', 'S${subdiv}_n${numIntervals}_i${interval}_low.png']},\n";
        }
        print FLIPPER "];";
        close(FLIPPER);
        print "    ['$numIntervals Intervals:Interval $interval', '$flipperPath'],\n";
    }
}
print "];\n";

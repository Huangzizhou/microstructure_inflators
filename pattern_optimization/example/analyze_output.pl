#!/usr/bin/perl -w
use strict;

my (@targetModuli, @currentModuli, @moduliDiff, @JS, @gradPNorm);
while (<>) {
    if (/Target moduli \(tensor\):\t(.*)/) {
        @targetModuli = split("\t", $1);
    }
    if (/Current moduli:\t(.*)/) {
        push(@currentModuli, [split("\t", $1)]);
    }
    if (/J_S = (.*)/) {
        push(@JS, $1);
    }
    if (/\|\|grad_p\|\|:\t(.*)/) {
        push(@gradPNorm, $1);
    }
}

my $maxIter = $#currentModuli;
# print(join("\t", @targetModuli) . "\n");

for my $i (0..$maxIter) {
    print($i . "\t"); # 1
    for my $m (0..$#targetModuli) { print($currentModuli[$i][$m] . "\t"); } # 2-5
    for my $m (0..$#targetModuli) { print(abs($targetModuli[$m] - $currentModuli[$i][$m]) . "\t"); } # 6-9
    print($JS[$i] . "\t"); # 10
    print($gradPNorm[$i]); # 11
    for my $m (0..$#targetModuli) { print("\t" . $targetModuli[$m]); } # 12-15
    print("\n");
}

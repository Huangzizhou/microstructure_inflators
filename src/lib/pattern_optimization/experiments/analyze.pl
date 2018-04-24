#!/usr/bin/env perl
use strict;
use warnings;

((@ARGV == 2) || (@ARGV == 3)) || die("Usage: analyze.pl pattern_num output_prefix [dir]\n");
my $pattern = $ARGV[0];
my $outputPrefix = $ARGV[1];

open(PRINTABLE_OUT, "> $outputPrefix.printable.txt") or die($!);
open(UNPRINTABLE_OUT, "> $outputPrefix.unprintable.txt") or die($!);

if (@ARGV == 2) {
    chdir $ARGV[1];
}

# open(OUTPUT, ">&STDOUT");
while (<exp_*.txt>) {
    open(RESULT_FILE, "< $_") or die("Couldn't open $_\n");
    my $numParams;
    my $iterationParsed = 0;
    my ($E, $nu, $anisotropy, $printable, @params);
    while (<RESULT_FILE>) {
        if (/^p:\t(.*)/) {
            @params = split(/\s+/, $1); 
            chomp(@params);
            $numParams //= @params;
            if ($numParams != @params) { die("Parameter number mismatch"); }
        }
        if (/^moduli:\t(.*)$/) {
            my @moduli = split(/\s+/, $1);
            chomp(@moduli);
            $E = $moduli[0];
            $nu = $moduli[3];
        }
        if (/^anisotropy:\t(.*)/) {
            $anisotropy = $1;
            chomp($anisotropy);
        }
        if (/^printable:\t(.*)/) {
            $printable = $1;
            chomp($printable);
            $iterationParsed = 1;
        }
        if ($iterationParsed) {
            if ($printable == 1) {
                print PRINTABLE_OUT "$pattern\t$E\t$nu\t$anisotropy\t$numParams\t" . join("\t", @params) . "\n";
            }
            else {
                print UNPRINTABLE_OUT "$pattern\t$E\t$nu\t$anisotropy\t$numParams\t" . join("\t", @params) . "\n";
            }
            $iterationParsed = 0;
        }
    }
    close(RESULT_FILE);
}
close(PRINTABLE_OUT);
close(UNPRINTABLE_OUT);

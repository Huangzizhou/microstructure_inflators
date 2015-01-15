#!/usr/bin/env perl
# Converts the dofs in a lookup table to the full set of dof values (without any
# equality constraints).
use strict;
use warnings;

(@ARGV == 1) || die("Usage: convert_lut_to_unconstrained_dof.pl lut.txt\n");

my $lutFile = $ARGV[0];
open(LUT, '<', $lutFile);

my %equalityConstraints;
open(EQCONST, '<', "equality_constraints.txt");
while (<EQCONST>) {
    my @row = split(/\t/);
    chomp(@row);
    my $pattern = $row[0];
    my @slice = @row[1..$#row];
    $equalityConstraints{$pattern} = \@slice;
}

while (<LUT>) {
    # pattern young poisson anisotropy #dofs dof0  ...
    my @row = split(/\t/);
    if ((@row < 5) || (@row != 5 + $row[4])) { die("Invalid row size\n"); }
    my $pattern = $row[0];
    my @constraint;
    if (exists $equalityConstraints{$pattern}) {
        @constraint = @{$equalityConstraints{$pattern}};
    }
    else { die("Couldn't look up pattern $pattern in equality_constraints.txt\n"); }
    my @params = @row[5..$#row];
    my $paramString = join("\t", @params);
    my $micro = ${ENV{'MICRO_DIR'}};
    my $constraintString = "";
    for my $c (@constraint) {
        $constraintString .= " -C'$c'"
    }
    my $fullParamString = `$micro/pattern_optimization/Inflator_cli -I -F$constraintString -p'$paramString' $micro/patterns/3D/reference_wires/pattern$pattern.wire | grep 'Full params'`;
    chomp($fullParamString);
    $fullParamString =~ s/Full params:\t//;
    my @fullParams = split(/\t/, $fullParamString);
    my $numFullParams = @fullParams;

    print(join("\t", @row[0..3]) . "\t$numFullParams\t$fullParamString\n");
}

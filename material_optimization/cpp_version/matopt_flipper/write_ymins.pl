#!/usr/bin/env perl
# Examine all runs matching a given pattern, determining a common ymin value for
# both objective and gradient norm plots, and writing these values to .ymin file
# for each run.
use strict; use warnings;
use File::Path qw(make_path);
use File::Find;
use List::Util qw(min max);

(@ARGV > 2) && die("Usage: write_ymins.pl [directory [filter]]\n");

my $rootDir = (@ARGV > 0) ? $ARGV[0] : '.';
my $filter  = (@ARGV > 1) ? $ARGV[1] : '.*';

my @matchingMSH;
find(sub { my $name = $File::Find::name;
           if ($name =~ m#$filter# && $name =~ m#\.msh$#) {
               push(@matchingMSH,$name) }
     }, $rootDir);
my (@allEnergies, @allGradNorms);
for my $file (@matchingMSH) {
    (my  $dataFile = $file) =~ s/\.msh$/.txt/;
    open(my $run, '<', $dataFile) or die("Couldn't open $dataFile for reading\n");
    while (<$run>) {
        if (/^([0-9]+) objective, gradient norm:\t(\S*)\t(\S*)$/) {
            push(@allEnergies, $2);
            push(@allGradNorms, $3);
        }
    }
}

my $yminString = join(' ', (min(@allEnergies), min(@allGradNorms)));

for my $file (@matchingMSH) {
    (my $yminFile = $file) =~ s/\.msh$/.ymin/;
    open(my $ymf, '>', $yminFile) or die ("couldn't open $yminFile for writing\n");
    print $ymf "$yminString\n";
    close($ymf);
    print "$yminFile\n";
}

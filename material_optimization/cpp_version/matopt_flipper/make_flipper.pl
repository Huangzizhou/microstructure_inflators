#!/usr/bin/perl -w
# Generate a flipper for a .msh file generated by the cpp material optimization
# code.
use strict;
use Cwd 'abs_path';
use File::Basename;
my $scriptDir = abs_path(dirname(__FILE__));

(@ARGV != 3 && @ARGV != 4) && die("Usage: make_flipper.pl matopt_out.txt fields.msh problem_name [matprops]\n");
my $optOutputFile = $ARGV[0];
my $mshFile = $ARGV[1];
my $problemName = $ARGV[2];
my $matPropString = (@ARGV == 4) ? $ARGV[3] : "E nu";

open(my $output, "<", $optOutputFile);
my (@energies, @gradNorms);
while (<$output>) {
    if (/^([0-9]+) objective, gradient norm:\t(\S*)\t(\S*)$/) {
        ($1 != scalar @energies) && die("Iterations must be 0-based, sequential\n");
        push(@energies, $2);
        push(@gradNorms, $3);
    }
}
open(my $pltData, ">", 'tmp.dat');
my $lastIt = $#energies;
for my $i (0..$lastIt) { print $pltData ("$i\t${energies[$i]}\t${gradNorms[$i]}\n"); }

################################################################################
# Render Fields
################################################################################
my $imagePrefix = $problemName;
my @fieldNames = ("Optimized u", "Neumann u", "Dirichlet u", "Young&apos;s Modulus", "Poisson Ratio");
my @fields = ("u", "u_neumann", "u_dirichletTargets", split(' ', $matPropString));
my $drawCalls = join("\\\n", map(qq(field="$_"; Call DrawField;), @fields));
`cat $scriptDir/render.geo | sed 's/<NITER>/$lastIt/; s/<PREFIX>/$imagePrefix/; s/<DRAW_CALLS>/$drawCalls/' > tmp.geo`;
`gmsh -n $mshFile tmp.geo 2>/dev/null`;
unlink 'tmp.geo';

################################################################################
# Trim whitespace
################################################################################
# for my $field (@fields) {
#     for my $i (0..$lastIt) {
#         my $image = "$imagePrefix.$i.$field.png";
#         `convert -trim $image $image`;
#     }
# }


# Gnuplot
my @titles = ("Objective", "Gradient Norm");
my @names =  ("objective", "gradNorm");
for my $col (2..2 + $#names) {
    for my $i (0..$lastIt) {
        my $plotImage = "$imagePrefix.$i.${names[$col-2]}.png";
        open(my $GP, "| gnuplot") or die("Couldn't open pipe to GnuPlot: $!\n");
        print $GP <<GNU_EOF;
        set term pngcairo size 800,600;
        set output '$plotImage';
        set xlabel 'Iteration';
        set logscale y;
        set yrange [1e-13:1];
        set autoscale ymax;
        set title '${titles[$col-2]}';
        unset key;
        set grid;
        plot 'tmp.dat' using 1:$col with lines, \\
             '< cat tmp.dat | grep -e "^$i\\s"' using 1:$col with points lc rgb 'black' ps 2 pt 6 notitle;
GNU_EOF
        close($GP);
    }
}

################################################################################
# Generate frames file.
################################################################################
my @views=(@fieldNames, "Objective Value", "Gradient Norm");
my @imageNames=(@fields, "objective", "gradNorm");
my @statistics=("objective", "gradient Norm");
my @statisticsData=(\@energies, \@gradNorms);

my $viewList = join(', ', map(qq('$_'), @views));
my $statList = join(', ', map(qq('$_'), @statistics));

print <<HERE;
title = 'Material optimization &quot;$problemName&quot;';
statistics = ['name', $statList];
views = [$viewList];
frames = [
HERE

for my $i (1..$lastIt) {
    print "\t{ name: 'Iter $i', 'image': [";
    print join(', ', map(qq!'$imagePrefix.$i.$_.png'!, @imageNames));
    print "]";
    for my $s (0..$#statistics) {
        my $val = ${$statisticsData[$s]}[$i];
        print ", '${statistics[$s]}': '$val'";
    }
    print(($i < $lastIt) ? "},\n" : "}\n");
}

print "];\n";

unlink 'tmp.dat';

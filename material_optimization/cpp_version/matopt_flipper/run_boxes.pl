#!/usr/bin/env perl
# Sweeps over James' 2D and 3D box examples, running material optimization.
# TODO: add filters for mesh/bc name
use strict;
use Time::HiRes;

(@ARGV < 1) && die("Usage: run_boxes.pl resultsDirectory [dim]\n");
my $resultsDir = $ARGV[0];
my @dims = (@ARGV == 2) ? ($ARGV[1]) : (2, 3);

my @regWeights3D = qw(0.0 0.01 0.1 1.0 10.0);
my @regWeights2D = qw(0.0 0.00001 0.0001 0.001 0.01 0.1 1.0);
my @regWeights = (\@regWeights2D, \@regWeights3D);

mkdir $resultsDir;

my $totalTime = 0;
for my $dim (@dims) {
    for my $meshPath (glob "box_compression/${dim}D/meshes/*") {
        my $mesh = $meshPath;
        $mesh =~ s#^.*?([^/]+)\.[^.]+$#\1#;
        for my $bcPath (glob "box_compression/${dim}D/BC/*/mesh_specific/$mesh.bc") {
            my $cond = $bcPath;
            $cond =~ s#.*/BC/([^/]+)/mesh_specific/.*#\1#;
            for my $material qw(orthotropic isotropic) {
                for my $regWeight (@{$regWeights[$dim - 2]}) {
                    my $name = "${dim}D:$mesh:$cond:$material:$regWeight";
                    my $start = Time::HiRes::time();
                    `$ENV{'MATOPT_DIR'}/MaterialOptimization_cli $meshPath $bcPath $resultsDir/$name.msh --numIters 20 --regularizationWeight=$regWeight --material=$material > $resultsDir/$name.txt`;
                    my $elapsed = Time::HiRes::time() - $start;
                    $totalTime += $elapsed;
                    print("$name\t$elapsed\n");
                }
            }
        }
    }
}

print "Total time:\t$totalTime\n";

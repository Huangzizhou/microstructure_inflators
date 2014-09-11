#!/usr/bin/env perl
# Sweeps over James' 2D and 3D box examples, running material optimization.
# TODO: add filters for mesh/bc name
use strict;
use Time::HiRes;
use File::Path qw(make_path);

(@ARGV < 1) && die("Usage: run_boxes.pl resultsDirectory [dim]\n");
my $resultsDir = $ARGV[0];
my @dims = (@ARGV == 2) ? ($ARGV[1]) : (2, 3);

my @itersPerDirichlet = qw(1 2 4 8 16 32);
my @anisoPenalty = qw(0.0 0.001 0.01);
my @regWeights3D = qw(0.0 0.01 0.1);
my @regWeights2D = qw(0.0 0.000001 0.00001 0.0001 0.001 0.01);
my @regWeights = (\@regWeights2D, \@regWeights3D);

my $totalTime = 0;

sub optimize {
    my ($dim, $meshPath, $mesh, $bcPath, $regWeight, $itersPerDirichlet, $material, $anisoPenalty) = @_;

    my $cond = $bcPath;
    $cond =~ s#.*/BC/([^/]+)/mesh_specific/.*#\1#;

    my $outpath = "$resultsDir/${dim}D/$mesh/$cond/$material/IPD_$itersPerDirichlet/AnisoPenalty_$anisoPenalty";

    make_path($outpath);
    $outpath .= "/reg_$regWeight";
    my $start = Time::HiRes::time();
    `$ENV{'MATOPT_DIR'}/MaterialOptimization_cli $meshPath $bcPath $outpath.msh --numIters 20 --regularizationWeight=$regWeight --material=$material --anisotropyPenaltyWeight $anisoPenalty > $outpath.txt`;
    my $elapsed = Time::HiRes::time() - $start;
    $totalTime += $elapsed;
    print("$outpath\t$elapsed\n");
}

for my $dim (@dims) {
    for my $meshPath (glob "box_compression/${dim}D/meshes/*") {
        my $mesh = $meshPath;
        $mesh =~ s#^.*?([^/]+)\.[^.]+$#\1#;
        for my $bcPath (glob "box_compression/${dim}D/BC/*/mesh_specific/$mesh.bc") {
            for my $regWeight (@{$regWeights[$dim - 2]}) {
                for my $itersPerDirichlet (@itersPerDirichlet) {
                    optimize($dim, $meshPath, $mesh, $bcPath, $regWeight, $itersPerDirichlet, 'isotropic', 0.0);
                    for my $anisoPenalty (@anisoPenalty) {
                        optimize($dim, $meshPath, $mesh, $bcPath, $regWeight, $itersPerDirichlet, 'orthotropic', $anisoPenalty);
                    }
                }
            }
        }
    }
}

print "Total time:\t$totalTime\n";

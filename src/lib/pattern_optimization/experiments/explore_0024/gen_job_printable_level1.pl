#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0024/printable/level1";
`mkdir -p $dir`;
my @inits = (
    "               0.3,                0.3,  -0.2204390029025167,    0.1974376094051841,                -0.3",
    "               0.3,                0.3,  -0.1898925750150332,  -0.00851096237812107, -0.2337904630697378",
    "0.3505113384324986,                0.3,                  0.3,    0.1486168333561043,                -0.3",
    "0.3544332145277191, 0.3410029150239251,  -0.1776739645992499, -0.004835228409135577, -0.2239153426911222",
    "0.3221060647120791, 0.3557906814458029, -0.03942524846503412,   -0.1867737573874208,                -0.3",
);
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3) {
        for my $radiusRange (0.6,0.8) {
            for my $poisson (0.20,0.225,0.25,0.275,0.30,0.325, 0.350,0.375,0.4,0.425) {
                for my $young (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_0024";
                    my $task = "exp_${young}_${poisson}_${radiusRange}_${transRange}_$i";
                    my $optFile = "$dir/$task.opt";
                    open(OPT_FILE, "> $optFile");
                    print OPT_FILE<<END;

{
	"dim": 3,
	"target": { "type": "isotropic",
                "young": $young,
                "poisson": $poisson },
    "initial_params": [${inits[$i]}],
	"radiusBounds": [0.3, $radiusRange],
	"translationBounds": [-$transRange, $transRange],
    "paramConstraints": ["p3 = p5 + 1.0/6.0"]
}
END
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern0024.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 4 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

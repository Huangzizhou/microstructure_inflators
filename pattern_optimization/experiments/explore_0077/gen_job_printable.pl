#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0077/printable";
`mkdir -p $dir`;
my @inits = ("0.3, 0.3,    0.1, -0.10",
             "0.3, 0.3,    0.1,  0.0",
             "0.3, 0.3,    0.15, 0.1",
             "0.3, 0.3,    0.2,  0.15",
             "0.3, 0.3,    0.2,  0.1");
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8,1.0) {
            for my $poisson (0.0,0.01,0.02,0.03,0.04,0.08,0.16,0.32) {
                for my $young (0.5, 1, 2, 4, 8, 16, 32, 64, 128) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_0077";
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
    "paramConstraints": ["p2 = p3"]
}
END
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern0077.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 16 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

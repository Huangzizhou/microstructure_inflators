#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0024/printable";
`mkdir -p $dir`;
my @inits = ("0.3, 0.3, 0.1, 0.0, -0.2",
             "0.5, 0.5, 0.1, 0.0, -0.2",
             "0.3, 0.3, -0.05, 0.0, -0.2"
);
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3) {
        for my $radiusRange (0.6,0.8) {
            for my $poisson (0.20,0.25,0.30,0.35,0.40,0.45) {
                for my $young (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8) {
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

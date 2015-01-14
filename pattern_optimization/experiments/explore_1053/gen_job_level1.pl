#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/1053/printable/level1";
`mkdir -p $dir`;
my @inits = (
    "               0.3, 0.3, 0.3, -0.03669997728327418,  0.2564440716127763, -0.04046528056313144",
    "0.3101236747231931, 0.3, 0.3, -0.05521141154434626,  0.0781276181503358,  -0.1007104258422012",
    "0.3336255663908697, 0.3, 0.3,                 -0.2, 0.04082332014080377,                 -0.2",
);

for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8) {
            for my $poisson (-0.15,-0.10,-0.05,0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.38) {
                for my $young (0.5, 1, 2, 3,  4, 8, 16, 32, 64) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_1053";
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
	"translationBounds": [-$transRange, $transRange]
}
END
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1053.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 4 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

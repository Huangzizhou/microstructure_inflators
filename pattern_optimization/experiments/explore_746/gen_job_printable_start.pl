#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0746/printable_fixed_gradient";
`mkdir -p $dir`;
my @inits = ("0.3, 0.3, 0.3, 0.025, -0.016666667, 0.025, -0.20");
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8,1.0) {
            for my $poisson (0.05,0.10,0.15,0.20,0.25,0.30,0.35) {
                for my $young (0.25, 1, 2, 4, 8, 16, 32, 64) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_746";
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
    "paramConstraints": ["p5 = p3 - 1/6"]
}
END
                    # my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1065.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material --dofOut $dir/${task}_params -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern0746.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 16 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

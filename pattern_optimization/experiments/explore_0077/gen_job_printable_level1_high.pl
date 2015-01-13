#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0077/printable/level1_high";
`mkdir -p $dir`;
my @inits = ("0.7581967362059817,                0.8,  0.1579464238313215,  0.07305407724734213",
             "0.6281393319728032,                0.8,   0.178653206216016,   0.1684896718257098",
             "0.5984510232723419, 0.8040160304872108,  0.1573005362600012, -0.02216170492765862");

for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.8,1.0) {
            for my $poisson (0.1,0.15,0.20,0.25,0.30) {
                for my $young (64, 128, 256) {
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
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 4 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

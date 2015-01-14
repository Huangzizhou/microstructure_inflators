#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0455/printable";
`mkdir -p $dir`;
# my @inits = (
#         "0.3, 0.3, 0.3,  0.1875, 0.02083400000000002,             0.104166, 0.1875, 0.02083400000000002, -0.02083400000000002",
#         "0.3, 0.3, 0.3,  0.0625, 0.02083400000000002, -0.02083400000000002, 0.1875, 0.02083400000000002,            -0.145834",
#         "0.3, 0.3, 0.3,  0.0625, 0.02083400000000002,             0.104166, 0.1875, 0.02083400000000002,            -0.145834",
#         "0.3, 0.3, 0.3, -0.0625, 0.02083400000000002, -0.02083400000000002, 0.1875, 0.02083400000000002, -0.02083400000000002",
#         "0.3, 0.3, 0.3, -0.0625, 0.02083400000000002,             0.104166, 0.1875, 0.02083400000000002, -0.02083400000000002",
#         "0.3, 0.3, 0.3, -0.0625,           -0.104166, -0.02083400000000002, 0.1875, 0.02083400000000002,  -0.02083400000000002"
# );
# The following are with parameter 6 eliminated using the equality constraint.
my @inits = (
        "0.3, 0.3, 0.3,  0.1875, 0.02083400000000002,             0.104166, 0.02083400000000002, -0.02083400000000002",
        "0.3, 0.3, 0.3,  0.0625, 0.02083400000000002, -0.02083400000000002, 0.02083400000000002,            -0.145834",
        "0.3, 0.3, 0.3,  0.0625, 0.02083400000000002,             0.104166, 0.02083400000000002,            -0.145834",
        "0.3, 0.3, 0.3, -0.0625, 0.02083400000000002, -0.02083400000000002, 0.02083400000000002, -0.02083400000000002",
        "0.3, 0.3, 0.3, -0.0625, 0.02083400000000002,             0.104166, 0.02083400000000002, -0.02083400000000002",
        "0.3, 0.3, 0.3, -0.0625,           -0.104166, -0.02083400000000002, 0.02083400000000002,  -0.02083400000000002"
);
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8,1.0) {
            for my $poisson (0.05,0.10,0.15,0.20,0.25,0.30,0.35, 0.40, 0.45) {
                for my $young (0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_0455";
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
    "paramConstraints": ["p6 = p7 - 1.0/6.0"]
}
END
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern0455.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 4 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}


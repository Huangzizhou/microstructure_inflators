#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0455/printable/level1";
`mkdir -p $dir`;

my @inits = (
"               0.3,                0.3,                0.3,                 -0.2, -0.03806803966507945,   0.1428049709483793,                   -0.2,                 0.2",
"               0.3,                0.3,                0.3,  -0.1170149735147282, 0.009982269677520382, 0.004292845537769796, -0.0003084098515694252, 0.09607697427948446",
"0.4243159072712671, 0.4082510842437387, 0.5463591274804483,   0.1826643799960756,   0.1213160595778835,  0.01660227885405416,     0.0957344675044378, 0.01881968820943329",
"0.3293052385189282, 0.4108764266275968,                0.6,  0.03771789699495515,  -0.0380059003293313, -0.04542403836142352,   -0.07403193429468598,  0.1123614086545771",
"0.3265757063029321, 0.7416025694592346, 0.7950515052866249, -0.03043997752871211, -0.05814743346568072,  0.01756848700299685,     -0.062199895719581,  0.1466178512448684",
" 0.433126333007766, 0.4075110146588668, 0.5980317029985058,  -0.0074020397847006, -0.04204925504316866,  -0.1041968466890823,   -0.06513244243940507, 0.09043945301285751"
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


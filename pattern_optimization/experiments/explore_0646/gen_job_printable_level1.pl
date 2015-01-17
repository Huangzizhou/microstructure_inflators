#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0646/printable/level1";
`mkdir -p $dir`;
my @inits = (
"               0.3, 0.3052459298491861,                0.3, -0.06686396055660668, 0.08873890519410477,                -0.4",
"               0.3,                0.3,                0.3, -0.05265998350940355, 0.06001180647958501,                -0.4",
"0.4315269691845966, 0.3466525141621321, 0.3522855833828589, -0.04805941386854557, 0.05331969358587606, -0.3264984113798033",
"0.5617740685756919, 0.4218640113654619, 0.3624406698572351, -0.01593325715795444,  0.1460360230948598, -0.2711623029285558"
);

    
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4,0.5) {
        for my $radiusRange (0.5, 0.6) {
            for my $poisson (0.10,0.15, 0.175, 0.20,0.25,0.30) {
                for my $young (0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_0646";
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
    "paramConstraints": ["p3 = p5 + 1.0/6.0",
                         "p4 = p7 + 1.0/6.0"]
}
END
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern0646.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 4 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

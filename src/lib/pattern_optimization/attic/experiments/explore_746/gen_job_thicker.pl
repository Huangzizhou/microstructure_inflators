#!/usr/bin/env perl
# Generates job with a higher thickness upper bound for DoF 2 (printability reasons).
# Also uses the B9Creator material.
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0746/printable_thicker_dof2_B9";
`mkdir -p $dir`;
my @inits = (
    "0.3848094605931864,                0.3, 0.4714045596016497, -0.03189375685172282,                  0.2, -0.1473645416166646,                -0.2",
    " 0.452813052133847,                0.3, 0.5350770233910694, 0.001271732267333516,                  0.2, -0.1447132359452263,                -0.2",
    "0.3000105650617761,                0.3, 0.6778426139309857,  0.04345112628608877,  0.06236888872928548,  0.2265049091791909, -0.1513064809098471",
    "               0.3,                0.3, 0.5531971034283635,  0.02325797885892794, -0.04993334050177621, 0.06537884598163084, -0.1511765298192638",
    " 0.402861244302206, 0.3355691672836764, 0.4863997996521016,  0.06907108113608314,  0.02602872811004467,  0.0337205845401634, -0.1549184169601092",
    "0.5476925224940168, 0.3529659145315275, 0.5999591847701783,  0.09091621975829335,  0.03003453868974218, 0.05371754676345092, -0.1728546991622326"
);

for (my $i = 0; $i < @inits; $i++) {
    for my $lowerDoFBound (0.375, 0.40, 0.45) {
        for my $transRange (0.2,0.3,0.4) {
            for my $radiusRange (0.6,0.8,1.0) {
                for my $poisson (0.2) {
                    for my $young (0.1, 0.2, 0.4, 0.8, 1.6) {
                        my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_746";
                        my $task = "exp_${young}_${poisson}_${radiusRange}_${transRange}_${lowerDoFBound}_$i";
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
    "paramConstraints": ["p5 = p3 - 1/6"],
    "bounds": { "var": 2, "lower": $lowerDoFBound, "upper": $radiusRange }
}
END
                        my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern0746.wire -m ~/microstructures/materials/B9Creator.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                        `create_pbs_noredirect.sh "$task" "$cmd" 2 4 1 0 > $dir/$task.pbs`;
                    }
                }
            }
        }
    }
}

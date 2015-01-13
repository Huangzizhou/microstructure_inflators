#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0077/printable/level1";
`mkdir -p $dir`;
my @inits = (
    "               0.3,                0.3,   0.144523348782698,   0.1316735460348962",
    "               0.3, 0.3026338002428404,  0.1445882522766346,  0.08237376920289716",
    "               0.3, 0.5497661975433412, 0.03597244727437442,  0.02453414325948702",
    "0.3047222542483074,                0.8, 0.04627996583313067,  0.02460450094218469",
    "0.3108653759733649,                  1, 0.09929762449612413,   0.0465317727627316",
    "0.3094250990348407, 0.9841192897880527,  0.2178936563095829,  0.09785819578270773",
    " 0.519269768703541,                0.3,  0.1690653711234929,  0.05306079275699779"
);

for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8,1.0) {
            for my $poisson (0.0,0.01,0.02,0.04,0.08,0.15,0.20,0.25, 0.30) {
                for my $young (0.5, 1, 2, 4, 8, 16) {
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

#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/1053/printable/level2";
`mkdir -p $dir`;
my @inits = (
"              0.3,                0.3,                0.3,                 -0.2,  0.0405631663574642,                -0.2",
"              0.3,                0.3,                0.3,   -0.191574860212959, 0.02729652628670853, -0.2402481203628092",
"              0.3,                0.3,                0.3,  -0.1002516053095339,  0.1216628150098014, -0.1026136732495584",
"0.305423339834501, 0.3039204072599164, 0.3543421283643307,  -0.0960032890426626,  0.1096249265711775, -0.1111946258565579",
"              0.3,                0.3, 0.3052692475165682, -0.05987310950362501,                 0.2, -0.0634921967051077"
);

for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8) {
            for my $poisson (-0.225, -0.20,-0.18, -0.15,-0.10,-0.05,0.0,0.05,0.1) {
                for my $young (2, 2.5,, 3, 3.5,, 3.75,  4, 8, 16) {
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

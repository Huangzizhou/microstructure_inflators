#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso";
`mkdir -p $dir`;

my @inits = ("0.3365684447812214,     0.3399461517585384,     0.85              ,     0.1820600831951494  ,   -0.07299983975961456,   0.3",
             "0.3051913175874095,     0.530919660465258 ,     0.3267816702045935,     -0.05093804176672936,   -0.2351652107224184 ,   0.296027017069701",
             "0.6393186248001879,     0.5988038483703614,     0.3438092691350659,     -0.06501055758437294,   -0.193016471920528  ,   0.2",
             "0.4904829985894038,     0.6053535473510607,     0.4902738000015927,     0.06074287842672824 ,   -0.1786010674550154 ,   0.3",
             "0.3               ,     0.6490949066679608,     0.6307104379252843,     0.2335993435548414  ,   -0.1525219643024219 ,   0.2946514469994678",
             "0.6740922495605617,     0.5558465957465085,     0.5377602947083551,     -0.07776134453164839,   -0.1284765263354063 ,   0.3",
             "0.8170553702994046,     0.7078306918266096,     0.849944802644253 ,     0.1834992967020583  ,   -0.1011521129709589 ,   0.2514000921255241"
             );
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $poisson (0.29,0.28,0.27,0.26) {
            for my $thickness (0.7,0.85,0.9) {
                for my $young (0.75,1..60) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_1065";
                    my $task = "exp_${young}_${poisson}_${transRange}_${thickness}_$i";
                    my $optFile = "$dir/$task.opt";
                    open(OPT_FILE, "> $optFile");
                    print OPT_FILE<<END;

{
	"dim": 3,
	"target": { "type": "isotropic",
                "young": $young,
                "poisson": $poisson },
    "initial_params": [${inits[$i]}],
	"radiusBounds": [0.3, $thickness],
	"translationBounds": [-$transRange, $transRange]
}
END
                    # my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1065.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material --dofOut $dir/${task}_params -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1065.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs.sh "$task" "$cmd" 2 16 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

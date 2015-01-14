#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso";
`mkdir -p $dir`;
my @inits = (" 0.5, 0.5, 0.5, 0.125, -0.125, 0.125",
             "0.6048051028377729, 0.8265562682840638, 0.5846618953987952, 0.2, -0.2, 0.0144785936511406",
             "0.3081734742450845, 0.7806828263155368,  0.6964524096341878,     -0.003007556001880967,     -0.2077176988443565, -0.0298784721100029",
             "0.3402013493117534, 0.6572802982771672,  0.85,                                     0.2,     -0.1619957275200183, 0.2",
             "0.3011123388001799, 0.5577363819200681,  0.8019609235918423,       0.02635485116573926,     -0.1843027719350475, 0.2",
             "0.3,                0.752282932439747,   0.6428481750918853,       0.3,                     -0.1845316509763536, 0.1648345183196891");
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3) {
        for my $poisson (0.3,0.34,0.40) {
            for my $young (1..60,65,70,75,80) {
                my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_1065";
                my $task = "exp_${young}_${poisson}_${transRange}.$i";
                my $optFile = "$dir/$task.opt";
                open(OPT_FILE, "> $optFile");
                print OPT_FILE<<END;

{
	"dim": 3,
	"target": { "type": "isotropic",
                "young": $young,
                "poisson": $poisson },
    "initial_params": [${inits[$i]}],
	"radiusBounds": [0.3, 0.85],
	"translationBounds": [-$transRange, $transRange]
}
END
                # my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1065.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material --dofOut $dir/${task}_params -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1065.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material --I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                `create_pbs.sh "$task" "$cmd" 2 16 1 0 > $dir/$task.pbs`;
            }
        }
    }
}

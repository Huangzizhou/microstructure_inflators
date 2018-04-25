#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso";
`mkdir -p $dir`;

my @inits = ("0.3004952743749066, 0.3233831440479957, 0.85, 0.1968913565046608, -0.125524701238178, 0.2",
             "0.4400845524524327, 0.5585012803944822, 0.4241066314500382, 0.2999999550642772, -0.1471482203620302, 0.3",
             "0.3435260202290866, 0.6603295901566544, 0.8411046559450476, 0.2, -0.1634192817626721, 0.2",
             "0.3201765666163743, 0.6463847577883233, 0.8413993895181772, 0.02326297032746096, -0.1865661749481731, 0.2",
             "0.3722214930500226, 0.6965433121789274, 0.8359451933148127, 0.09636983153856915, -0.1770025510704293, 0.2",
             "0.7372307978836318, 0.85, 0.8495620645049942, 0.1999989144981233, -0.1088645949255852, 0.09450278538120285");
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4,0.5) {
        for my $poisson (0.32,0.30,0.25,0.20) {
            for my $young (0.75,1..60,65,70,75,80,100) {
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
                my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1065.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                `create_pbs.sh "$task" "$cmd" 2 16 1 0 > $dir/$task.pbs`;
            }
        }
    }
}

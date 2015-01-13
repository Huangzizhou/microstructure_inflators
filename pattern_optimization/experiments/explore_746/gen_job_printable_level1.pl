#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0746/printable_level1";
`mkdir -p $dir`;
my @inits = ("0.3448461104082935, 0.3, 0.337063857658326, 0.0467532169813901, 0.1379563763332891, -0.07033184704312299, -0.1195409224436434",
             "0.3332845771630324,                0.3, 0.3944463197820134,  0.1288284265410754,    0.106962165183772,  -0.0128771119675714, -0.0881306867415543",
             "0.5695762076148453,                0.3, 0.6679283676136766,  0.1224742564343336,   0.1279026181918073, -0.04131894950191675, -0.1032685944691676",
             "0.3530760211975006, 0.3001093722158493, 0.4062902272570716, 0.02343073535181012,  0.03913585799394298,   0.0120303436714766, -0.1541009690012572",
             "0.4393276044346218, 0.3033072557341753, 0.5121018331236895, 0.04839651437778981,  0.02266146717135321,  0.03266519195922183, -0.1781889963400446",
             "0.5586351925479098, 0.3645703679314365, 0.5571797460358554,  0.0602182730343838, -0.02028599119675316,  0.08387177316506775,                -0.2",
             "0.3110016754941763, 0.3153733070523612, 0.3580749764795764, 0.04929333906613562,  0.01697159600133206,  0.02874971649679563, -0.1582522342617083",
             "0.3065298427823533, 0.3025078887449232, 0.3212235041383399, 0.04847128540886213, -0.01212356939235678,  0.03648505797589275, -0.1922608249211767");
for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8,1.0) {
            for my $poisson (0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40) {
                for my $young (0.25, 1, 2, 4, 8, 16, 32, 64) {
                    my $exploreDir = "~/microstructures/pattern_optimization/experiments/explore_746";
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
    "paramConstraints": ["p5 = p3 - 1/6"]
}
END
                    # my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern1065.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material --dofOut $dir/${task}_params -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    my $cmd = "~/microstructures/pattern_optimization/PatternOptimization_cli $optFile -p $exploreDir/pattern0746.wire -m ~/MeshFEM/experiments/fit_validation/ProJet7000_2D.material -I -S1 -A simple -v0 --solver levenberg_marquardt > $dir/$task.txt 2>&1";
                    `create_pbs_noredirect.sh "$task" "$cmd" 2 4 1 0 > $dir/$task.pbs`;
                }
            }
        }
    }
}

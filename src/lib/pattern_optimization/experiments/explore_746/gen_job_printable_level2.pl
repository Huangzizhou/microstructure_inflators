#!/usr/bin/env perl
use strict;
my $dir ="${ENV{'SCRATCH'}}/explore_iso/0746/printable/level2";
`mkdir -p $dir`;
my @inits = ("0.3, 0.3, 0.3, -0.02389871202846582, 0.2223613136455214, -0.1604724418532668, -0.3",
             "0.3038853251029391,                0.3, 0.3457399235182927,  0.1341928942441754,   0.1330199550292794, -0.02012803670763464,  -0.0323043663999593",
             " 0.313873019183756, 0.3169189188226065,                0.3,  0.1137616035192222,   0.1326080176027487, -0.02842219428480425, -0.06085131359456381",
             "0.5241777619910103, 0.3165079194743983, 0.6921009817854828,  0.0927508473008037,  0.01527116632588905,  0.02369645383119874, -0.08860589209521323",
             "0.3518646101827675,  0.375325379867986,                0.3, 0.03384382044007544, -0.03848911955260712,  0.03161009545402195,  -0.2384009233220629"
);

for (my $i = 0; $i < @inits; $i++) {
    for my $transRange (0.2,0.3,0.4) {
        for my $radiusRange (0.6,0.8,1.0) {
            for my $poisson (0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40) {
                for my $young (0.25, 0.5, 0.75, 1, 2, 4, 8) {
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

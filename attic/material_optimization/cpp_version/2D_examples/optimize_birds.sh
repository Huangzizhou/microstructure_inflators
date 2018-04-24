# usage: ./optimize_birds.sh [resultsDir]
resultsDir=${1-results}
dir=$resultsDir/bird_opt
mkdir -p $dir

# $1: mesh, $2: bc, $3: name, $4: args
function run {
    echo "$MeshFEM/MaterialOptimization_cli $resultsDir/$1 $2 $4 $dir/$3.msh | tee $dir/$3.txt" > $dir/$3.sh
    sh $dir/$3.sh
}

run  bird_tri.msh          bcs/flap_up.bc          narrow:flap_up:bird      "-b var_bounds/narrow.js -R -r0.00000001 -n50"
run bird2_tri.msh          bcs/flap_up.bc          narrow:flap_up:bird2     "-b var_bounds/narrow.js -R -r0.000001   -n50"
run  bird_tri.msh          bcs/flap_mixed.bc       narrow:flapmix:bird      "-b var_bounds/narrow.js -R -r0.00000001 -n50"
run  bird_wingless_tri.msh bcs/flap_up_wingless.bc narrow:flap_up:wingless  "-b var_bounds/narrow.js    -r0.0000001  -n16"
run bird_wingless2_tri.msh bcs/flap_up_wingless.bc narrow:flap_up:wingless2 "-b var_bounds/narrow.js    -r0.0000001  -n16"
run bird_stubby_tri.msh    bcs/flap_up_stubby.bc   narrow:flap_up:stubby    "-b var_bounds/narrow.js    -r0.000001   -n30 -N3"

run  bird_tri.msh          bcs/flap_up.bc          moderate:flap_up:bird      "-b var_bounds/moderate.js -R -r0.00000001 -n50"
run bird2_tri.msh          bcs/flap_up.bc          moderate:flap_up:bird2     "-b var_bounds/moderate.js -R -r0.000001   -n50"
run  bird_tri.msh          bcs/flap_mixed.bc       moderate:flapmix:bird      "-b var_bounds/moderate.js -R -r0.00000001 -n50"
run  bird_wingless_tri.msh bcs/flap_up_wingless.bc moderate:flap_up:wingless  "-b var_bounds/moderate.js    -r0.0000001  -n16"
run bird_wingless2_tri.msh bcs/flap_up_wingless.bc moderate:flap_up:wingless2 "-b var_bounds/moderate.js    -r0.0000001  -n16"
run bird_stubby_tri.msh    bcs/flap_up_stubby.bc   moderate:flap_up:stubby    "-b var_bounds/moderate.js    -r0.000001   -n30 -N3"

run  bird_tri.msh          bcs/flap_up.bc          wide:flap_up:bird      "-b var_bounds/wide.js -R -r0.00000001 -n50"
run bird2_tri.msh          bcs/flap_up.bc          wide:flap_up:bird2     "-b var_bounds/wide.js -R -r0.000001   -n50"
run  bird_tri.msh          bcs/flap_mixed.bc       wide:flapmix:bird      "-b var_bounds/wide.js -R -r0.00000001 -n50"
run  bird_wingless_tri.msh bcs/flap_up_wingless.bc wide:flap_up:wingless  "-b var_bounds/wide.js    -r0.0000001  -n16"
run bird_wingless2_tri.msh bcs/flap_up_wingless.bc wide:flap_up:wingless2 "-b var_bounds/wide.js    -r0.0000001  -n16"
run bird_stubby_tri.msh    bcs/flap_up_stubby.bc   wide:flap_up:stubby    "-b var_bounds/wide.js    -r0.000001   -n30 -N3"

# usage: ./optimize_bars.sh [resultsDir]
resultsDir=${1-results}
dir=$resultsDir

# $1: mesh, $2: bc, $3: name, $4: args
function run {
    echo "$MeshFEM/MaterialOptimization_cli $resultsDir/$1 $2 $4 $dir/$3.msh | tee $dir/$3.txt" > $dir/$3.sh
    sh $dir/$3.sh
}

for r in {0..2}; do
    mkdir $dir/buckle_smoother_refine_$r;
    for i in {half_period,half_period.75N,bump_half,bump_third,flat,full_period}; do
        run bar_2D_4x20.$r.msh bcs/bar_2D.$i.bc buckle_smoother_refine_$r/$i "-b var_bounds/bar.js -r 0.000001 -n25"
    done;
done

for r in {0..2}; do
    mkdir $dir/buckle_sharper_refine_$r;
    for i in {half_period,half_period.75N,bump_half,bump_third,flat,full_period}; do
        run bar_2D_4x20.$r.msh bcs/bar_2D.$i.bc buckle_sharper_refine_$r/$i "-b var_bounds/bar.js -r 0.0000001 -n25"
    done;
done

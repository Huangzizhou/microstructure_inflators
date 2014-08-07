#!/usr/bin/env zsh
# Example of how to run a sweep over parameters of the material optimization
# code on various inputs.
# Assumes the MATOPT_DIR variable holds the name of the cpp material
# optimization directory
[[ -d $1 ]] || { echo "Usage: matopt_run_example.sh example_dir [results_dir]"; exit -1; }
mesh_dir=$1;
total=0.0;
dir=${2-results};
mkdir $dir;
for regWeight in 0.0 0.00001 0.0001 0.001 0.01 0.1 1.0; do
    for mesh in $mesh_dir/*.obj; do
        mesh_name=$(basename $mesh .obj)
        for cond in $mesh_dir/$mesh_name*.bc; do
            cond_name=${$(basename $cond .bc)#$mesh_name.}
            out_name=${mesh_name}:${cond_name}:reg$regWeight
            startTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
            $MATOPT_DIR/TestMaterialOptimization2D $mesh $cond 20 $regWeight > $dir/$out_name.txt
            mv mtest.msh $dir/$out_name.msh;
            endTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
            elapsed=$((($endTime - $startTime) / 1000.0));
            total=$(($total + $elapsed))
            echo "$out_name:\t" $(printf "%0.3fs" $elapsed);
        done;
    done
done
echo "Total time:\t$total";

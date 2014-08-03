#!/usr/bin/env zsh
# Example of how to run a sweep over parameters of the material optimization
# code on various inputs.
# Assumes the MATOPT_DIR variable holds the name of the cpp material
# optimization directory
[[ -d $1 ]] || { echo "Usage: matopt_run_example.sh example_dir"; exit -1; }
mesh_dir=$1;
total=0.0;
for regWeight in {10,100,1000,10000}; do
    dir=mtest_$regWeight
    mkdir $dir;
    for mesh in $mesh_dir/*.obj; do
        for cond in ${mesh%.obj}*.bc; do
            name="${cond%bc}"
            startTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
            $MATOPT_DIR/TestMaterialOptimization2D $mesh $cond 15 $regWeight > $dir/$(basename $name).txt
            mv mtest.msh $dir/$(basename $name).msh;
            endTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
            elapsed=$((($endTime - $startTime) / 1000.0));
            total=$(($total + $elapsed))
            echo "$name:\t" $(printf "%0.3fs" $elapsed);
        done;
    done
done
echo "Total time:\t$total";

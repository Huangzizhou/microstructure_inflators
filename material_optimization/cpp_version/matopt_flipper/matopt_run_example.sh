# Example of how to run a sweep over parameters of the material optimization
# code on various inputs.
total=0.0;
for regWeight in {10,100,1000,10000}; do
    dir=mtest_$regWeight
    mkdir $dir;
    cd $dir;
    for size in {8x8,16x16,32x32}; do
        for cond in {tilt,compress,full_period,half_period}; do
            startTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
            ../TestMaterialOptimization2D ~/Downloads/mesh_bc/square_$size.obj ../box_bc/${cond}_$size.bc 15 $regWeight > ${cond}_$size.txt;
            mv mtest.msh ${cond}_$size.msh;
            endTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
            elapsed=$((($endTime - $startTime) / 1000.0));
            total=$(($total + $elapsed))
            echo "${cond}_$size:\t" $(printf "%0.3fs" $elapsed);
        done;
    done
    cd ..;
done
echo "Total time:\t$total";

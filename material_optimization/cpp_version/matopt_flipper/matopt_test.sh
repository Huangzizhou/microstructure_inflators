total=0.0;
for bounds in {bounded,unbounded}; do
    for regWeight in {10,100,1000,10000}; do
        dir=mtest_${bounds}_$regWeight
        mkdir $dir;
        cd $dir;
        for size in {8x8,16x16,32x32}; do
            for cond in {tilt,compress,full_period,half_period}; do
                startTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
                ../TestMaterialOptimization2D_$bounds ~/Downloads/mesh_bc/square_$size.obj ../box_bc/${cond}_$size.bc 15 $regWeight > ${cond}_$size.txt;
                mv mtest.msh ${cond}_$size.msh;
                endTime=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`
                elapsed=$((($endTime - $startTime) / 1000.0));
                total=$(($total + $elapsed))
                echo "${cond}_$size:\t" $(printf "%0.3fs" $elapsed);
            done;
        done
        cd ..;
    done
done
echo "Total time:\t$total";

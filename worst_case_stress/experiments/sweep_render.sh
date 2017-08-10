result_dir=$1
for run in $(lfs find $result_dir -name '*.msh' | sed 's/it_//' s'/_[0-9][0-9]*.remeshed.msh//' | sort -u); do
    mkdir -p $run
    echo "{'cmd': '/home/fjp234/microstructures/worst_case_stress/experiments/make_iterate_animation.sh $run'},"
done
for msh in $(lfs find $result_dir -name '*.msh'); do
    name=$(basename $msh .remeshed.msh)
    rundir=$(echo $msh | sed 's/it_//' s'/_[0-9][0-9]*.remeshed.msh//')
    num=${name##*_}
    echo "gmsh_offscreen $msh $MICRO_DIR/worst_case_stress/view_ptwise_measure.opt -o $rundir/$num.png"
done

dir=$(timestamped_path.sh $SCRATCH/pattern_optimization/random/brick5)
./make_jobs.sh brick5 brick5_base.opt $MICRO_DIR/patterns/3D/brick5.wire $MeshFEM/experiments/fit_validation/ProJet7000_2D.material simple levenberg_marquardt $dir 2 32 16 0
echo $dir

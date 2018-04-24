dir=$(timestamped_path.sh $SCRATCH/gradient_validation/brick5)
./make_jobs.sh brick5 ../random_opt/brick5_base.opt $MICRO_DIR/patterns/3D/brick5.wire $MeshFEM/experiments/fit_validation/ProJet7000_2D.material simple 7 -0.02 0.02 $dir 2 32 24 0
./make_jobs.sh brick5 ../random_opt/brick5_base.opt $MICRO_DIR/patterns/3D/brick5.wire $MeshFEM/experiments/fit_validation/ProJet7000_2D.material simple 7  -0.2  0.2 $dir 2 32 24 0
./make_jobs.sh brick5 ../random_opt/brick5_base.opt $MICRO_DIR/patterns/3D/brick5.wire $MeshFEM/experiments/fit_validation/ProJet7000_2D.material loop 7 -0.02 0.02 $dir 2 32 24 0
./make_jobs.sh brick5 ../random_opt/brick5_base.opt $MICRO_DIR/patterns/3D/brick5.wire $MeshFEM/experiments/fit_validation/ProJet7000_2D.material loop 7  -0.2  0.2 $dir 2 32 24 0
echo $dir

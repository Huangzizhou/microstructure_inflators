$MICRO_DIR/worst_case_stress/WCSOptimization_cli -p $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj -m $MICRO_DIR/materials/B9Creator.material \
    target_tensor_job.opt  -M 2d_meshing_opts.opt --ortho_cell --vertexThickness \
    --WCSWeight 1e-300 --JSWeight 1.0 --TensorFitConstraint \
    --solver slsqp -o it

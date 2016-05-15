#!/usr/bin/env zsh
mkdir -p fine_adaptive_meshes coarse_meshes
for smooth in {15,20,30,60}; do
    $MICRO_DIR/worst_case_stress/GradientComponentValidation -R -P10 --JSWeight 1024 -m $MICRO_DIR/materials/B9Creator.material -p $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj <(./octacell_smooth_job.sh $smooth)  -IV -M <($MICRO_DIR/worst_case_stress/gradient_validation/octacell_iso/mesh_options_adaptive.sh 1e-03 500) -d1 0 --range_relative 0.0 1 -s 0 -o out > /dev/null
    mv out_0.msh fine_adaptive_meshes/smooth_$smooth.msh

    $MICRO_DIR/worst_case_stress/GradientComponentValidation -R -P10 --JSWeight 1024 -m $MICRO_DIR/materials/B9Creator.material -p $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj <(./octacell_smooth_job.sh $smooth)  -IV -M <($MICRO_DIR/worst_case_stress/gradient_validation/octacell_iso/mesh_options_adaptive.sh 1e-03 192) -d1 0 --range_relative 0.0 1 -s 0 -o out > /dev/null
    mv out_0.msh coarse_meshes/smooth_$smooth.msh
done

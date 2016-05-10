#!/usr/bin/env zsh
mkdir -p fine_meshes
mkdir -p coarse_meshes
for smooth in {robust,weak}; do
    $MICRO_DIR/worst_case_stress/GradientComponentValidation -m $MICRO_DIR/materials/B9Creator.material -p $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj $smooth.opt -IV -M <(./mesh_options_adaptive.sh 1e-03 4096) -d1 0 --range_relative 0.0 1 -s 0 -o out > /dev/null
    mv out_0.msh fine_meshes/$smooth.msh

    $MICRO_DIR/worst_case_stress/GradientComponentValidation -m $MICRO_DIR/materials/B9Creator.material -p $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj $smooth.opt -IV -M <(./mesh_options_adaptive_forced_bdrylen.sh 5e-01 384 0.06) -d1 0 --range_relative 0.0 1 -s 0 -o out
    # $MICRO_DIR/worst_case_stress/GradientComponentValidation -m $MICRO_DIR/materials/B9Creator.material -p $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj $smooth.opt -IV -M <(./mesh_options_adaptive.sh 1e-01 384) -d1 0 --range_relative 0.0 1 -s 0 -o out
    mv out_0.msh coarse_meshes/$smooth.msh
done

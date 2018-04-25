#!/usr/bin/env zsh
 $MICRO_DIR/worst_case_stress/WCSOptimization_cli --JSWeight 0.0 --WCSWeight 1 \
    -p $MICRO_DIR/patterns/3D/reference_wires/pattern1065.wire \
    -m $MICRO_DIR/materials/B9Creator.material -V --solver slsqp \
    $MICRO_DIR/worst_case_stress/experiments/resolution_test/young_0.35.opt \
    -M coarser_3d_meshing_opts.opt \
    -O -C -P8 -R -n200


# TODO: have this call make_iterate_animation and do the plotting.
run=$1
gnuplot -e "run='$run'" $MICRO_DIR/worst_case_stress/bp_convergence_plots.gpi
$MICRO_DIR/worst_case_stress/experiments/make_iterate_animation.sh $run

result_dir=$1
for run in $(lfs find $result_dir -name '*.msh' | sed 's/it_//; s/_[0-9.][0-9.]*.remeshed.msh//' | sort -u); do
    gnuplot -e "run='$run'" $MICRO_DIR/worst_case_stress/bp_convergence_plots.gpi
done

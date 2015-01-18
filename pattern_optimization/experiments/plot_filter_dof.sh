#!/usr/bin/env zsh
if [[ ! $# -eq 2 ]]; then
    echo "usage: plot_filter_dof.sh data.txt num_thickness_dofs"
    exit;
fi
data=$1;
numThicknessDofs=$2
for dof in {0..$((numThicknessDofs - 1))}; do
    # DoF's field in the lookup table field--awk fields are 1-indexed!
    dofField=$(($dof + 5 + 1));
    for lowerBound in {0.3,0.325,0.35,0.375,0.40,0.425,0.45,0.475}; do
gnuplot <<HERE;
set term png size 800,600;
unset key;
set grid;
set ylabel "Young's Modulus";
set xlabel "Poisson Ratio";
set yrange [1.0/100:200];
set xrange [-0.30:0.5];
set logscale y;
set output "plot_${dof}_${lowerBound}.png";
set title "Lower Bound $lowerBound on DoF $dof";
plot "< $HOME/microstructures/pattern_optimization/experiments/filter_field.sh $data $dofField $lowerBound 9999999" using 3:2 with points lc rgb "red" pt 1;
HERE
    done
done

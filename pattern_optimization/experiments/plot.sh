# run after analyze.sh
if [[ ! $# -eq 3 ]]; then
    echo "usage: plot.sh isotropy_dist printable.txt unprintable.txt"
    exit;
fi
dist=$1;
printable=$2
unprintable=$3
gnuplot <<HERE;
set term png size 1024,768;
unset key;
set grid;
set ylabel "Young's Modulus";
set xlabel "Poisson Ratio";
set yrange [1.0/100:200];
set xrange [-0.30:0.5];
set logscale y;
set output "plot_all_$dist.png";
set title "Anisotropy in Range 1 +/- $dist, All";
plot "< $HOME/microstructures/pattern_optimization/experiments/filter.sh $dist $unprintable" using 3:2 with points lc rgb "red" pt 2, \
     "< $HOME/microstructures/pattern_optimization/experiments/filter.sh $dist $printable" using 3:2 with points lc rgb "green" pt 1;
set title "Anisotropy in Range 1 +/- $dist, Printable Only";
set output "plot_printable_$dist.png";
plot "< $HOME/microstructures/pattern_optimization/experiments/filter.sh $dist $printable" using 3:2 with points lc rgb "black" pt 1;
HERE

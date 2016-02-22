if [ $# -ne 1 ]; then
    echo "usage: plot_all.sh results_dir"
    exit;
fi;
resultDir=$1
for p in {1..6}; do for param in {0..13}; do gnuplot -e "resultDir='$resultDir'; P=$p; paramNum=$param" plot.gpi; done; done

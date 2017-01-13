SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cp $SCRIPT_DIR/{analyze.py,data_at_stress_reduction.py,plot.gpi} .
mkdir filtered_data/
for thresh in {0.005,0.01,0.02,0.05,0.1}; do
    python analyze.py $thresh $thresh > filtered_data/filtered_$thresh.txt;
done

mkdir plots
gnuplot plot.gpi
mv coverage*.png plots

pushd plots
cp $SCRIPT_DIR/frames.js .
~/Research/flipper/install.sh
popd

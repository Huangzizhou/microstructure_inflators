for deg in {1,2}; do for type in {direct,projected}; do for p in {0..8}; do gnuplot -e "nsvType='$type'; deg =$deg; paramNum=$p" plot.gpi; done; done; done

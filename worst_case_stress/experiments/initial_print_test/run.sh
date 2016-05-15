# ./gen_meshes.sh
# ./tile.sh
./gen_data_traction.sh | tee traction_experiment.dat
./gen_data_dirichlet.sh | tee dirichlet_experiment.dat
# gnuplot plot.gpi
./extrude.sh
# ./render.sh

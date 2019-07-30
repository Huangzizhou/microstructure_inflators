#!/bin/bash

BIN_FOLDER="/Users/davitozoni/Desktop/NYU/Research/Repositories/microstructures/cmake-build-release/"

for f in png/*.png; 
do 
  name=`echo $f | sed 's/png\///' | sed 's/\.png//'`
  echo $name

  ../../../cmake-build-release/voxel2Mesh --mopts 2d_meshing_opts.json $f "msh/"$name".msh"; 
  ../../../cmake-build-release/shape_optimization/polygonMesher_cli --periodic "msh/"$name".msh" "msh/"$name"_coarse.msh" -m 2d_coarse_meshing_opts.json

  mkdir "initial_"$name

  $BIN_FOLDER/worst_case_stress/WCSOptimization_cli -i BoundaryPerturbation -p "msh/"$name"_coarse.msh" -m Default.material target_tensor_job.opt --vertexThickness --WCSWeight 1.0 --proximityRegularizationWeight 1.0 --smoothingRegularizationWeight 10.0 --solver custom_gradient_descent -o "initial_"$name"/it" -n 1 --usePthRoot --pnorm 12
  cp "initial_"$name"/it_sol.txt" "results_msh/"$name"_initial.txt"
  cp "initial_"$name"/it_best" "results_msh/"$name"_initial.msh"
  
  rm -rf "initial_"$name
done
 

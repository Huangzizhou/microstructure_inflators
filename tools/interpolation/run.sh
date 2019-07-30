#!/bin/bash

BIN_FOLDER="/Users/davitozoni/Desktop/NYU/Research/Repositories/microstructures/cmake-build-release/"

# restart
#rm -r msh
#rm -r results_msh
#rm -r results_png

mkdir msh
mkdir results_msh
mkdir results_png

# copy initial and target file
first_file=`ls png/ | sort -n | head -1`
first_name=`echo $first_file | sed 's/png\///' | sed 's/\.png//'`
last_file=`ls png/ | sort -n | tail -1`
last_name=`echo $last_file | sed 's/png\///' | sed 's/\.png//'`
cp "png/"$first_file "results_png/"$first_name"..png"
cp "png/"$last_file "results_png/"$last_name"f.png"

# also, compute volume of initial and target structures
../../../cmake-build-release/voxel2Mesh --mopts 2d_meshing_opts.json png/$first_file "msh/"$first_name"...msh"; 
../../../cmake-build-release/voxel2Mesh --mopts 2d_meshing_opts.json png/$last_file "msh/"$last_name"f.msh"; 

target1_volume=`../../../../MeshFEM/cmake-build-release/MeshFEM/PeriodicHomogenization_cli "msh/"$first_name"...msh" --outputVolume | grep Volume | head -1 | awk '{print $2}'`
target2_volume=`../../../../MeshFEM/cmake-build-release/MeshFEM/PeriodicHomogenization_cli "msh/"$last_name"f.msh" --outputVolume | grep Volume | head -1 | awk '{print $2}'`

# compute number of frames
number_frames=`ls -l all_png/ | grep -v tile | tail -n +2 | wc -l`

# loop through all pngs in folder pngs
#  1- transform png into msh (and put it into folder msh)
#  2- optimize msh to reduce stress (and put result into folder msh_result)
#  3- transform msh into png (and put it into folder png_result)
id=0
for f in png/*.png; 
do
  if [[ $f == *"tiled"* ]]; then
    continue
  fi
  
  name=`echo $f | sed 's/png\///' | sed 's/\.png//'`
  echo $name
  initial_mesh="msh/"$name".msh"

  if [ -f "results_msh/"$name".msh" ]
  then
    continue
  fi
  
  ../../../cmake-build-release/voxel2Mesh --mopts 2d_meshing_opts.json $f "msh/"$name".msh"; 
  ../../../cmake-build-release/shape_optimization/polygonMesher_cli --periodic $initial_mesh "msh/"$name"_coarse.msh" -m 2d_coarse_meshing_opts.json
  cp "msh/"$name"_coarse.msh" "msh/"$name".msh"

  t=`python -c "print(1.0*$id/($number_frames-1))"`
  initial_volume=`python -c "print($t * $target1_volume + (1.0-$t) * $target2_volume)"`
  #initial_volume=`../../../../MeshFEM/cmake-build-release/MeshFEM/PeriodicHomogenization_cli "msh/"$name".msh" --outputVolume | grep Volume | head -1 | awk '{print $2}'`
  initial_volume="1.464"
  #echo "t: "$t
  #echo "initial_volume: "$initial_volume
  for ((it=1;it<=10;it++))
  do
    iterations_folder="iterations_"$name"_"$it
    mkdir $iterations_folder

    #../../../cmake-build-release/shape_optimization/polygonMesher_cli --periodic $initial_mesh "msh/"$name"_coarse.msh" -m 2d_coarse_meshing_opts.json --collapseCleanup
    ../../../cmake-build-release/shape_optimization/polygonMesher_cli --periodic $initial_mesh "msh/"$name"_coarse.msh" -m 2d_coarse_meshing_opts.json

    ../../../cmake-build-release/worst_case_stress/GenBoundaryPerturbationInflatorJob -o 0.1 --targetVolume $initial_volume "msh/"$name"_coarse.msh" > job.json

    $BIN_FOLDER/worst_case_stress/WCSOptimization_cli -i BoundaryPerturbation -p "msh/"$name"_coarse.msh" -m Default.material job.json --targetVolWeight 1000.0 --vertexThickness --WCSWeight 1.0 --WCSTarget 20.0 --proximityRegularizationWeight 0.5 --smoothingRegularizationWeight 10000.0 --solver custom_bfgs -o "$iterations_folder/it" -n 200 --usePthRoot --pnorm 12
    cp "$iterations_folder/it_best" "msh/"$name".msh"
  
  done

  cp "$iterations_folder/it_best" "results_msh/"$name".msh"
  cp "$iterations_folder/it_sol.txt" "results_msh/"$name".txt"
  ../../../../MeshFEM/cmake-build-release/MeshFEM/mesh_convert -b "results_msh/"$name".msh" "results_msh/"$name".off"
  ./off2png "results_msh/"$name".off" "results_png/"$name".png"
  
  rm -rf "iterations_"$name"_"*

  id=`python -c "print($id + 1)"`
done
 

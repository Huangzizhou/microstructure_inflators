#!/bin/bash

for f in results_msh/*_initial.msh; 
do 
  name=`echo $f | sed 's/results_msh\///' | sed 's/\_initial.msh//'`
  optimized_file="results_msh/"$name".msh"

  ../../../../MeshFEM/cmake-build-release/MeshFEM/PeriodicHomogenization_cli $f -m Default.material --outputVolume > tmp_material.txt
  
  initial_mu=`cat tmp_material.txt | grep v_yx  | awk '{print $3}'`
  initial_E=`cat tmp_material.txt  | grep Young | awk '{print $4}'`
  initial_nu=`cat tmp_material.txt | grep shear | awk '{print $4}'`
  initial_vol=`cat tmp_material.txt | grep Volume | tail -1 | awk '{print $3}'`

  if [[ -f $optimized_file ]]
  then
    ../../../../MeshFEM/cmake-build-release/MeshFEM/PeriodicHomogenization_cli $optimized_file -m Default.material --outputVolume > tmp_material.txt
  
    mu=`cat tmp_material.txt | grep v_yx  | awk '{print $3}'`
    E=`cat tmp_material.txt  | grep Young | awk '{print $4}'`
    nu=`cat tmp_material.txt | grep shear | awk '{print $4}'`
    vol=`cat tmp_material.txt | grep Volume | tail -1 | awk '{print $3}'`

    echo $name $initial_mu $mu $initial_E $E $initial_nu $nu $initial_vol $vol
  fi
done



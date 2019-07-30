#!/bin/bash

for f in results_msh/*_initial.txt; 
do 
  name=`echo $f | sed 's/results_msh\///' | sed 's/\_initial.txt//'`
  optimized_file="results_msh/"$name".msh"

  if [[ -f $optimized_file ]]
  then

    initial_wcs=`cat $f | grep Ptw | sed 's/Max Ptwise WCS:	//'`
    result="results_msh/"$name".txt"
    wcs=`cat $result | grep Ptw | sed 's/Max Ptwise WCS:	//'`

    echo $name	$initial_wcs $wcs
  fi
done



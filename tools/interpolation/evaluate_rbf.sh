#!/bin/bash

n=$1
mkdir "rbf_"$n

for i in example_0001_0003/png/0001_to_0003_0[0-9][0-9].png
do 
  name=`echo $i | sed 's/\(.*\/\)\(.*\)\(\.png\)/\2/'`
  #error=`./png2levelset.py $i --basisPerDim $n "rbf_"$n"/"$name"_rbf.png" | grep Error | awk '{print $3}'`
  error_msg=`./png2levelset.py $i --basisPerDim $n "rbf_"$n"/"$name"_rbf.png" | grep "Error\|pixel"`
  #echo $error_msg
  error=`echo $error_msg | grep Error | awk '{print $3}'`
  wrong_pixels=`echo $error_msg | grep pixel | awk '{print $6}'`
  echo $name" "$error" "$wrong_pixels
done



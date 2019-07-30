#!/bin/bash

for f in png/*.png
do
  if [[ $f == *"tiled"* ]]; then
    continue
  fi

  name=`echo $f | sed 's/png\///' | sed 's/\.png//'`
  echo $name
  ./tile_pngs.sh $f "png/"$name"_tiled.png" 3 3 1 1
done

for f in results_png/*.png
do
  if [[ $f == *"tiled"* ]]; then
    continue
  fi

  name=`echo $f | sed 's/results_png\///' | sed 's/\.png//'`
  echo $name
  ./tile_pngs.sh $f "results_png/"$name"_tiled.png" 3 3 1 1
done

convert -delay 20 -loop 0 png/*_tiled.png initial.gif
convert -delay 20 -loop 0 results_png/*_tiled.png optimized.gif



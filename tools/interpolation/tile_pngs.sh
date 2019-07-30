#!/bin/bash

image="$1"
result="$2"
num_rows="$3"
num_cols="$4"
focus_cell_x="$5"
focus_cell_y="$6"


# generate focused cell
convert $image PNG8:focused.png
convert focused.png -depth 12 -fuzz 10% -fill 'rgb(255,0,0)' -opaque black focused.png

# construct focused column
cp $image $result
for ((i=1;i<$num_rows;i++))
do
  if [ $i == $focus_cell_y ]
  then
    convert $result focused.png -append $result
  else
    convert $result $image -append $result
  fi
done
cp $result tmp_focused_rows.png

cp $image $result
# construct other columns
for ((i=1;i<$num_rows;i++))
do
  convert $result $image -append $result
done
cp $result tmp_rows.png

for ((j=1;j<$num_cols;j++))
do
  if [ $j == $focus_cell_x ]
  then
    convert $result tmp_focused_rows.png +append $result
  else
    convert $result tmp_rows.png +append $result
  fi
done

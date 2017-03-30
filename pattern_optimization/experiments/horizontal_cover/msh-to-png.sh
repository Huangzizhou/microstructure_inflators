#!/bin/bash

echo "Print \"$2\"; Exit;" > print-script.txt

gmsh $1 print-script.txt



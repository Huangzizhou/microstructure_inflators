#!/usr/bin/env zsh
mesh=$1
python cellfaceCoordAnalysis.py <(./coord_extract $mesh) 0 1.0 | head

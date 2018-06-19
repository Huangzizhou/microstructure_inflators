#!/bin/bash
../ShapeOptimization_cli -b dirichlet.json -p octa_cell_square.obj -m ../../materials/B9Creator.material sample_job.opt -M 2d_meshing_opts.opt --solver gradient_descent --pnorm 10 --usePthRoot --nIters 50 --step 0.001 -o it

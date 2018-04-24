#!/bin/bash
../../ShapeOptimization_cli -b boundary_conditions.json -p inner_triangle_cantilever.wire -m ../../../materials/Russia.material inner_triangle_cantilever.job -M 2d_meshing_opts.opt --solver gradient_descent --pnorm 2 --usePthRoot --nIters 100 --step 0.001 -o iterations/it

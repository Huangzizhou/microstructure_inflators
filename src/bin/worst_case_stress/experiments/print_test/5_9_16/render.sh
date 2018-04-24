#!/usr/bin/env zsh
mkdir -p renders
for i in {15,20,30,60}; do
    gmsh -n smooth_${i}.msh render_mesh.opt
    mv render.png renders/smooth_${i}.mesh.png
    gmsh -n smooth_${i}.msh render_wcs.opt
    mv render.png renders/smooth_${i}.wcs.png
    gmsh -n smooth_${i}_pull_xy.0.msh render_pull_stress.opt
    mv render.png renders/smooth_${i}.pull_xy.png
    gmsh -n smooth_${i}_pull_y.0.msh render_pull_stress.opt
    mv render.png renders/smooth_${i}.pull_y.png
done

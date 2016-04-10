#!/usr/bin/env zsh
for cf in {coarse,fine}; do
    mkdir -p ${cf}_tilings
    for i in ${cf}_meshes/*.msh; do
        cp $i in.msh
        for r in {0..3}; do
            t=$((2**$r))
            cp in.msh ${cf}_tilings/$(basename $i .msh)_${t}x${t}.msh
            mesh_convert in.msh --Sx 0.5 --Sy 0.5 -r in.msh
        done
        rm in.msh
    done
done

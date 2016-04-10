#!/usr/bin/env zsh
mkdir -p sim_pull_y_traction sim_pull_xy_traction 
for mesh in fine_meshes/*.msh; do
    wcs=$(../../WCS_cli $mesh -m $MICRO_DIR/materials/B9Creator.material -d2 | cut -f2)
    echo -ne "$i\t$wcs\t"
    for tmesh in fine_tilings/$(basename $mesh .msh)_{1,2,4,8}*.msh; do
        pully_mesh=sim_pull_y_traction/$(basename $tmesh)
        pullxy_mesh=sim_pull_xy_traction/$(basename $tmesh)
        
        Simulate_cli -d2 $tmesh -b pull_y.bc  -m $MICRO_DIR/materials/B9Creator.material -o $pully_mesh > /dev/null
        Simulate_cli -d2 $tmesh -b pull_xy.bc -m $MICRO_DIR/materials/B9Creator.material -o $pullxy_mesh > /dev/null

        pully_stress=$(msh_processor  $pully_mesh  -e 'stress' -l --maxMag --max)
        pullxy_stress=$(msh_processor $pullxy_mesh -e 'stress' -l --maxMag --max)

        echo -ne "$pully_stress\t$pullxy_stress\t"
    done
    echo ""
done

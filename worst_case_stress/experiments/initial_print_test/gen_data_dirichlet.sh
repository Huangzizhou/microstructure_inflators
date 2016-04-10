#!/usr/bin/env zsh
mkdir -p sim_pull_y_dirichlet sim_pull_xy_dirichlet 
for mesh in fine_meshes/*.msh; do
    i=${$(basename $mesh .msh)##smooth_}
    wcs=$(../../WCS_cli $mesh -m $MICRO_DIR/materials/B9Creator.material -d2 | cut -f2)
    echo -ne "$i\t$wcs\t"
    for tmesh in fine_tilings/$(basename $mesh .msh)_{1,2,4,8}*.msh; do
        # Stretch all tilings with a unit force.
        yForce=$(Simulate_cli -d2 $tmesh -b pull_y_dirichlet.bc  -m $MICRO_DIR/materials/B9Creator.material -o dummy.msh | grep "region 2" | cut -f3;)
        yDispMag=$((1.0 / $yForce))
        xyForce=$(Simulate_cli -d2 $tmesh -b pull_xy_dirichlet.bc  -m $MICRO_DIR/materials/B9Creator.material -o dummy.msh | grep "region 2" | cut -f3;)
        xyDispMag=$((1.0 / $xyForce))
        rm dummy.msh

        pully_mesh=sim_pull_y_dirichlet/$(basename $tmesh)
        pullxy_mesh=sim_pull_xy_dirichlet/$(basename $tmesh)
        
        Simulate_cli -d2 $tmesh -b <(./pull_y_dirichlet_scaled.sh $yDispMag)   -m $MICRO_DIR/materials/B9Creator.material -o $pully_mesh > /dev/null
        Simulate_cli -d2 $tmesh -b <(./pull_xy_dirichlet_scaled.sh $xyDispMag) -m $MICRO_DIR/materials/B9Creator.material -o $pullxy_mesh > /dev/null

        pully_stress=$(msh_processor  $pully_mesh  -e 'stress' -l --maxMag --max)
        pullxy_stress=$(msh_processor $pullxy_mesh -e 'stress' -l --maxMag --max)

        echo -ne "$pully_stress\t$pullxy_stress\t"
    done
    echo ""
done

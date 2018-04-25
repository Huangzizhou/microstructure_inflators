mkdir -p scaled_extrusions

# extrude coarse tilings
for i in extrusions/*.msh; do
    mesh_convert --Sx 10 --Sy 10 --Sz 200 $i scaled_extrusions/$(basename $i)
done

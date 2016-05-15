mkdir -p extrusions

# extrude coarse tilings
for i in coarse_tilings/*.msh; do
    mesh_convert $i -e 0.05 --quadAspectThreshold 80.0 extrusions/$(basename $i)
done

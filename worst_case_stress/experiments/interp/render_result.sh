cat <<HEADER_END
title = 'Pattern interpolation experiment';
views = ['render'];
statistics = ['moduli', 'anisotropy', 'params'];
framesLabel = 'step';
frames = [
HEADER_END

for i in step_{0..19}.txt; do
    imageName=${i%.txt}.png
    params=$(../../tools/analyze_constrained_stdout.py $i ProximityRegularization 8e-9 | grep params | cut -f2-)
    moduli=$(../../tools/analyze_constrained_stdout.py $i ProximityRegularization 8e-9 | grep moduli | cut -f2-)
    anisotropy=$(../../tools/analyze_constrained_stdout.py $i ProximityRegularization 8e-9 | grep anisotropy | cut -f2-)
    echo "{'image': ['$imageName'], 'moduli': '$moduli', 'anisotropy': '$anisotropy', 'params': '$params'},"
    # ../../../isosurface_inflator/isosurface_cli 2D_orthotropic $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj -p "$params" -m <(../../gradient_validation/octacell_iso/mesh_options.sh 1e-4 2048) tmp.msh > /dev/null
    # gmsh_offscreen tmp.msh view.opt -o tmp.png 2> /dev/null > /dev/null
    # ~/animations/colorize_gmsh.sh tmp.png $imageName
done
echo "];"
# rm tmp.png
# trim_image_sequence.py step_*.png

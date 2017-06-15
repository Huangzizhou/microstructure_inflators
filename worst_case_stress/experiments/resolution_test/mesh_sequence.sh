script_dir=$MICRO_DIR/worst_case_stress/experiments/resolution_test/

# Test on 1065
# Mesh at the optimum to decide on the sequence of mesh resolutions
params="0.6073229928107067      0.7242407847863863      0.217254184891536 0.7536703422371146      0.7532888389419082      0.7537447253219028 0.620307344960584       0.2233776811593922    0.7253559147271538 0.6178134056489273      0.2259006162733315      0.7172171372016858 0.07565026600749267     0.04681935440518284     0.05296816488690493 0.04072742792297104     0.05099841848242969    0.05052533462443665 0.04745099274255049     0.0546643729278777      0.02254334455374395 0.05939162207292249     0.069435576998302       0.06261029204416511 0.05490183776660545    0.07189059276715995     0.05513635900514388 0.07084828789376167"

for i in {0..6}; do
    $MICRO_DIR/isosurface_inflator/isosurface_cli orthotropic $MICRO_DIR/patterns/3D/reference_wires/pattern1065.wire -m <( $script_dir/mesh_opts.py $i) \
        --params "$params" -O mopts_$i.msh
done
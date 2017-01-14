pat=$1
num=$2
dir=$3

# dir is, e.g. $SCRATCH/wcsmin/autocover_coarser_p8_1_12_17

script_dir=$MICRO_DIR/worst_case_stress/experiments/min_autocover

pattern=$MICRO_DIR/patterns/3D/reference_wires/pattern$pat.wire

optString=$(python $script_dir/extract_min_linf.py $dir/stdout_${pat}_$num.txt 0.05 0.05)
optParams=$(echo $optString | cut -f5-)

initParams=$(head -n $num $script_dir/tables/$pat.txt | tail -n1 | cut -f6-)

$MICRO_DIR/isosurface_inflator/isosurface_cli orthotropic $pattern \
    -m $script_dir/coarser_3d_meshing_opts.opt \
    --params "$optParams" $pat.$num.opt.msh

$MICRO_DIR/isosurface_inflator/isosurface_cli orthotropic $pattern \
    -m $script_dir/coarser_3d_meshing_opts.opt \
    --params "$initParams" $pat.$num.init.msh

optImg=${pat}_${num}.opt.png
initImg=${pat}_${num}.init.png
gmsh_offscreen $pat.$num.opt.msh  $script_dir/view_pat.opt -o $optImg
gmsh_offscreen $pat.$num.init.msh $script_dir/view_pat.opt -o $initImg

convert -trim -border 5x5 \
    $initImg $optImg +append \
    -gravity Center -background khaki \
    label:"Stress reduction, init stress, opt stress: $(echo $optString | cut -f2-4 | tr '\t' ',    ')" \
    -append ${pat}_$num.png

rm $pat.$num.opt.msh $pat.$num.init.msh $optImg $initImg

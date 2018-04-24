pat=$1
num=$2

mkdir -p ${pat}_${num}
cd ${pat}_${num}

script_dir=$MICRO_DIR/worst_case_stress/experiments/min_autocover/promising_for_experiments
pattern=$MICRO_DIR/patterns/3D/reference_wires/pattern$pat.wire

optParams=$(grep stdout_${pat}_$num.txt $script_dir/promising_db.txt | cut -f7-)
initParams=$(head -n $num $script_dir/../tables/$pat.txt | tail -n1 | cut -f6-)

optMesh=opt.msh
initMesh=init.msh

$MICRO_DIR/isosurface_inflator/isosurface_cli orthotropic $pattern \
    -m $script_dir/meshing_opts.opt -O \
    --params "$optParams" $optMesh

$MICRO_DIR/isosurface_inflator/isosurface_cli orthotropic $pattern \
    -m $script_dir/meshing_opts.opt -O \
    --params "$initParams" $initMesh

$MeshFEM/PeriodicHomogenization_cli -m $MICRO_DIR/materials/B9Creator.material -O $optMesh  > opt.homog.txt
$MeshFEM/PeriodicHomogenization_cli -m $MICRO_DIR/materials/B9Creator.material -O $initMesh > init.homog.txt

$MICRO_DIR/worst_case_stress/WCS_cli -m $MICRO_DIR/materials/B9Creator.material -O $optMesh  -o opt.wcs.msh
$MICRO_DIR/worst_case_stress/WCS_cli -m $MICRO_DIR/materials/B9Creator.material -O $initMesh -o init.wcs.msh

$MeshFEM/ConstStrainDisplacement_cli -m $MICRO_DIR/materials/B9Creator.material -O -s "0 0 -1 0 0 0" -f $optMesh  opt.zcompress.msh
$MeshFEM/ConstStrainDisplacement_cli -m $MICRO_DIR/materials/B9Creator.material -O -s "0 0 -1 0 0 0" -f $initMesh init.zcompress.msh

optCell=fullcell.$optMesh
initCell=fullcell.$initMesh
$MICRO_DIR/isosurface_inflator/replicate $optMesh $optCell
$MICRO_DIR/isosurface_inflator/replicate $initMesh $initCell
$MeshFEM/Simulate_cli -m $MICRO_DIR/materials/B9Creator.material $optCell  -b $script_dir/testing_setup.bc -o opt.zcompress.sim.msh
$MeshFEM/Simulate_cli -m $MICRO_DIR/materials/B9Creator.material $initCell -b $script_dir/testing_setup.bc -o init.zcompress.sim.msh

optTiling=tiling.$optMesh
initTiling=tiling.$initMesh
$MICRO_DIR/isosurface_inflator/replicate $optMesh  -r 2x2x1 $optTiling
$MICRO_DIR/isosurface_inflator/replicate $initMesh -r 2x2x1 $initTiling
$MeshFEM/Simulate_cli -m $MICRO_DIR/materials/B9Creator.material $optTiling  -b $script_dir/testing_setup.bc -o opt.zcompress.sim_2x2x1.msh
$MeshFEM/Simulate_cli -m $MICRO_DIR/materials/B9Creator.material $initTiling -b $script_dir/testing_setup.bc -o init.zcompress.sim_2x2x1.msh

rm $optMesh $initMesh $optCell $initCell $optTiling $initTiling

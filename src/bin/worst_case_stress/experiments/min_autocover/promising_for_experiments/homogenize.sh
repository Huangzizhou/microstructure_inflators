pat=$1
num=$2

cd ${pat}_${num}

$MeshFEM/PeriodicHomogenization_cli -m $MICRO_DIR/materials/B9Creator.material -O opt.wcs.msh  > opt.homog.txt
$MeshFEM/PeriodicHomogenization_cli -m $MICRO_DIR/materials/B9Creator.material -O init.wcs.msh > init.homog.txt

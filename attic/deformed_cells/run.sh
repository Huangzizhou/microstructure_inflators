# $1: mesh
# $2: param set name
# $2: parameters to try
mesh=`basename $1`;
name=$2
$MeshFEM/DeformedCells_cli $1 -m $MICRO_DIR/materials/B9Creator.material -p --homogenize < $3 > out.${mesh%.*}.$name.txt

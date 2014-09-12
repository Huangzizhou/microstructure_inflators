# Generate all the triangulated meshes from the source .obj geometry.
# usage: ./generate_geometry.sh [resultsDir]
resultsDir=${1-results}
tar -xf source_objs.tgz -C $resultsDir
cd $resultsDir

# Generate birds
$MeshFEM/mesh_convert bird.obj --Sx 5 --Sy 5 -q1 bird_tri.msh
$MeshFEM/mesh_convert bird.obj -s bird2.obj
$MeshFEM/mesh_convert bird2.obj --Sx 10 --Sy 10 -q1 bird2_tri.msh

$MeshFEM/mesh_convert bird_wingless.obj --Sx 5 --Sy 5 -q1 bird_wingless_tri.msh
$MeshFEM/mesh_convert bird_wingless.obj -s bird_wingless2.obj
$MeshFEM/mesh_convert bird_wingless2.obj --Sx 10 --Sy 10 -q1 bird_wingless2_tri.msh

$MeshFEM/mesh_convert bird_stubby.obj --Sx 5 --Sy 5 -q1 bird_stubby_tri.msh 

# Generate refined bars
for r in {0..2}; do
    $MeshFEM/mesh_convert bar_2D_4x20.obj -q$r bar_2D_4x20.$r.msh
done

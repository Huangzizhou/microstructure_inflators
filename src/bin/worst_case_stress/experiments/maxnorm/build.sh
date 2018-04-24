g++ -std=c++11 maxnorm.cc -I$CSGFEM -I$MeshFEM -I$EIGEN_INC -I$BOOST_INC $MeshFEM/MeshIO.cc $MeshFEM/Types.cc -o maxnorm

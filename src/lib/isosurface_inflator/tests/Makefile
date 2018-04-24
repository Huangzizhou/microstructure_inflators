all: coord_extract

coord_extract: coord_extract.cc
	g++ -std=c++11 coord_extract.cc $(MeshFEM)/MeshIO.cc -I$(MeshFEM) -I$(EIGEN_INC) -I${BOOST_INC} -o coord_extract

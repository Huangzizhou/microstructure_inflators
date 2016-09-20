# CFLAGS=-O0 -g
CFLAGS=-O2
all: data_viewer
data_viewer: data_viewer.cc $(MeshFEM)/MeshIO.o Makefile
	g++ $(CFLAGS) -std=c++11 -I$(MeshFEM) -I$(CSGFEM) -I$(EIGEN_INC) -I/opt/local/include -framework GLUT -framework OpenGL $^ -Wno-deprecated-declarations -L/opt/local/lib/ -lboost_filesystem-mt -lboost_system-mt -o data_viewer

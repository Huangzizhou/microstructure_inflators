all: LipschitzComponents
LipschitzComponents: LipschitzComponents.cc
	g++ -O2 -std=c++11 -I$(CSGFEM) -I$(EIGEN_INC) -I/opt/local/include -framework GLUT -framework OpenGL LipschitzComponents.cc -Wno-deprecated-declarations -L/opt/local/lib/ -lboost_filesystem-mt -lboost_system-mt -o LipschitzComponents

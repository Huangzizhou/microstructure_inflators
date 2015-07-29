all: data_viewer
data_viewer: data_viewer.cc
	g++ -O2 -std=c++11 -I/opt/local/include -framework GLUT -framework OpenGL data_viewer.cc -Wno-deprecated-declarations -L/opt/local/lib/ -lboost_filesystem-mt -lboost_system-mt -o data_viewer

CPPFLAGS=-O2
# CPPFLAGS=-O2 -DCGAL_MESH_3_VERBOSE
# CPPFLAGS=-g -O0 # -DCGAL_MESH_3_VERBOSE -DCGAL_MESH_3_PROTECTION_DEBUG -ftemplate-backtrace-limit=0
# CPPFLAGS += -DDEBUG_MARCHING_SQUARES
TARGETS=isosurface_cli
all: $(TARGETS)

isosurface_cli: isosurface_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o IsosurfaceInflator.o VCGSurfaceMesher.o BoxIntersectionMesher.o MeshingOptions.o CGALClippedVolumeMesher.o 
	g++ -L/usr/local/lib -L/opt/local/lib -lcgal -lboost_thread-mt -lboost_system-mt -lgmp -lmpfr -ladept -lboost_program_options-mt $^ -o $@
TestWireMesh: tests/TestWireMesh.o
	g++ -L/usr/local/lib $^ -o $@

INCLUDES+=-isystem /opt/local/include -I$(MeshFEM) -I$(CSGFEM) -isystem /opt/local/include/eigen3 -isystem $(VCGLIB_PATH) 
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DBENCHMARK
%.o: %.cc
	$(CXX) $(CPPFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(TARGETS)

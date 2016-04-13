MeshFEM ?= $(MESHFEM_PATH)
include $(MeshFEM)/platform_defs.mk

INCLUDES=-I$(MeshFEM) -I$(CSGFEM) -isystem $(EIGEN_INC) -isystem $(BOOST_INC) \
		 -isystem $(ADEPT_INC) -isystem $(VCGLIB_INC) \
		 -isystem $(CGAL_INC) -isystem $(GMP_INC) -isystem $(MPFR_INC) \
		 -isystem $(TRIANGLE_INC)
CPPFLAGS=-O2
# CPPFLAGS=-O0 -g

LIBS+=-L$(CGAL_LIB) -lCGAL -L$(GMP_LIB) -lgmp -L$(MPFR_LIB) -lmpfr
LIBS+=-L$(ADEPT_LIB) -ladept
LIBS+=$(BOOST_LFLAGS)
LIBS+=$(TRIANGLE_LFLAGS)

# CPPFLAGS=-O2 -DCGAL_MESH_3_VERBOSE
# CPPFLAGS=-g -O0 # -DCGAL_MESH_3_VERBOSE -DCGAL_MESH_3_PROTECTION_DEBUG -ftemplate-backtrace-limit=0
# CPPFLAGS += -DDEBUG_MARCHING_SQUARES
TARGETS=isosurface_cli SignedDistance_cli
all: $(TARGETS)

isosurface_cli: isosurface_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/GlobalBenchmark.o IsosurfaceInflator.o  MidplaneMesher.o MeshingOptions.o BoxIntersectionMesher.o # VCGSurfaceMesher.o CGALClippedVolumeMesher.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

SignedDistance_cli: SignedDistance_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

TestWireMesh: tests/TestWireMesh.o
	g++ -L/usr/local/lib $^ -o $@

INCLUDES+=-isystem /opt/local/include -I$(MeshFEM) -I$(CSGFEM) -isystem /opt/local/include/eigen3 -isystem $(VCGLIB_PATH) 
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DBENCHMARK
%.o: %.cc
	$(CXX) $(CPPFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(TARGETS)

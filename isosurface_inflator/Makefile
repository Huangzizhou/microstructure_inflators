MeshFEM ?= $(MESHFEM_PATH)
include $(MeshFEM)/platform_defs.mk

INCLUDES=-I$(MeshFEM) -I$(CSGFEM) -isystem $(EIGEN_INC) -isystem $(BOOST_INC) \
		 -isystem $(ADEPT_INC) -isystem $(VCGLIB_INC) \
		 -isystem $(CGAL_INC) -isystem $(GMP_INC) -isystem $(MPFR_INC) \
		 -isystem $(TRIANGLE_INC)

ifeq ($(DEBUG), 1)
	CPPFLAGS+=-O0 -g
else
	CPPFLAGS+=-O2
endif

ifeq ($(DNDEBUG), 1)
	CPPFLAGS+=-DNDEBUG
endif

LIBS+=-L$(CGAL_LIB) -lCGAL -L$(GMP_LIB) -lgmp -L$(MPFR_LIB) -lmpfr
LIBS+=-L$(ADEPT_LIB) -ladept
LIBS+=$(BOOST_LFLAGS)
LIBS+=$(TRIANGLE_LFLAGS)
LIBS+=$(TBB_LFLAGS)

# CPPFLAGS += -DCGAL_MESH_3_VERBOSE -DCGAL_MESH_3_DEBUG
# CPPFLAGS=-O2 -DCGAL_MESH_3_VERBOSE
# CPPFLAGS=-g -O0 # -DCGAL_MESH_3_VERBOSE -DCGAL_MESH_3_PROTECTION_DEBUG -ftemplate-backtrace-limit=0
# CPPFLAGS += -DDEBUG_MARCHING_SQUARES
TARGETS=isosurface_cli replicate PrintabilityCheck PrintabilityConstraints # SignedDistance_cli OffsetBoundsToTranslationBounds

# DEPENDENCY STUFF
DEPDIR:=.d
COMMA:=,
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS=-MT $@ -MMD -MP -MF $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).Td

all: $(TARGETS)

isosurface_cli: isosurface_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/GlobalBenchmark.o IsosurfaceInflator.o  MidplaneMesher.o MeshingOptions.o ConvexHullTriangulation.o BoxIntersectionMesher.o VCGSurfaceMesher.o CGALClippedVolumeMesher.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

replicate: replicate.o $(MeshFEM)/MeshIO.o $(MeshFEM)/MSHFieldParser.o $(MeshFEM)/Types.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

PrintabilityCheck: PrintabilityCheck.o $(MeshFEM)/MeshIO.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

PrintabilityConstraints: PrintabilityConstraints.o $(MeshFEM)/MeshIO.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

OffsetBoundsToTranslationBounds: OffsetBoundsToTranslationBounds.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/GlobalBenchmark.o IsosurfaceInflator.o  MidplaneMesher.o MeshingOptions.o BoxIntersectionMesher.o CGALClippedVolumeMesher.o VCGSurfaceMesher.o ../pattern_optimization/PatternOptimizationJob.o $(MeshFEM)/Materials.o ConvexHullTriangulation.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

SignedDistance_cli: SignedDistance_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o ConvexHullTriangulation.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

TestWireMesh: tests/TestWireMesh.o
	g++ -L/usr/local/lib $^ -o $@

INCLUDES+=-isystem /opt/local/include -I$(MeshFEM) -I$(CSGFEM) -isystem /opt/local/include/eigen3 -isystem $(VCGLIB_INC) 
INCLUDES+=$(TBB_IFLAGS)
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DBENCHMARK
CPPFLAGS+=$(WARNING_FLAGS)

# Delete built-in rule, forcing ours to run
%.o: %.cc
%.o: %.cc Makefile
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) -c $< -o $@
	mv -f $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).Td $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).d

%.o: %.c
%.o: %.c Makefile
	$(CC) $(DEPFLAGS) -c $(CFLAGS) $< -o $@
	mv -f $(DEPDIR)/$$(subst .,$(COMMA),$(subst /,__,$*)).Td $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).d

.PHONY: clean
clean:
	rm -f *.o $(TARGETS)

# Read dependencies
-include $(DEPDIR)/*.d

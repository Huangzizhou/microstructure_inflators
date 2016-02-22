include $(MeshFEM)/platform_defs.mk

INCLUDES=-I$(LEVMAR_INC) # bring levmar in first since we have a modified copy that might be overriden by the system's
INCLUDES+=-I$(MeshFEM) -I$(CSGFEM) -isystem $(EIGEN_INC) -I$(SUITESPARSE_INC) -I$(CLIPPER_INC) -I$(TRIANGLE_INC)  \
		  -isystem $(VCGLIB_INC) -isystem $(DLIB_INC)
INCLUDES+=-isystem $(CERES_INC) -isystem $(LIBMATHEVAL_INC)
INCLUDES+=-I$(PYMESH_PATH)/src -I$(PYMESH_PATH)/tools
INCLUDES+=-isystem ../Luigi/wireinflator2D/src

# stuff from patter_optimization
PATOPTIM=$(MICROSTRUCTURES_PATH)/pattern_optimization
INCLUDES+=-I$(PATOPTIM)

LIBS=$(BOOST_LFLAGS) $(SUITESPARSE_LFLAGS) $(CERES_LFLAGS) $(LIBMATHEVAL_LFLAGS)
LIBS+=$(CLIPPER_LFLAGS) -L$(PYMESH_PATH)/lib -lwires $(TRIANGLE_LFLAGS) 
LIBS+=$(LEVMAR_LFLAGS)


SOURCES=PatternOptimization_DeformedCell_cli.cc PatternOptimizationJob.cc\
		GetQuadDeformations_cli.cc  \
		TileQuadFromFile_cli.cc  FindCommonDeformations_cli.cc UnMap_cli.cc \
		GetQuadCenters_cli.cc ConvertToObj_cli.cc ExpandQuad_cli.cc testEdgeMesh.cc

DEFPATOPT_OBJS=PatternOptimization_DeformedCell_cli.o $(PATOPTIM)/PatternOptimizationJob.o \
            $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
            $(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o

GETQUADDEFS_OBJS=GetQuadDeformations_cli.o

MAP_OBJS=FindCommonDeformations_cli.o

UNMAP_OBJS=UnMap_cli.o

TILEQUAD_OBJS=TileQuadFromFile_cli.o \
            $(MeshFEM)/MeshIO.o \
			$(MeshFEM)/Types.o

# optional objects
testEdgeMesh_OBJS=testEdgeMesh.o \
            $(MeshFEM)/MeshIO.o \
			$(MeshFEM)/Types.o

EXPQUAD_OBJS=ExpandQuad_cli.o \
            $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o \
            $(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o

CONVERT_OBJS=ConvertToObj_cli.o \
            $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o

GETQUADCTRS_OBJS=GetQuadCenters_cli.o

OBJS=$(DEFPATOPT_OBJS) $(GETQUADDEFS_OBJS) $(TILEQUAD_OBJS) $(MAP_OBJS) $(UNMAP_OBJS) $(CONVERT_OBJS) $(EXPQUAD_OBJS) $(testEdgeMesh_OBJS) $(GETQUADCTRS_OBJS)

DEFAULT_TARGETS=PatternOptimization_DeformedCell_cli \
				GetQuadDeformations_cli \
				TileQuadFromFile_cli \
				FindCommonDeformations_cli UnMap_cli
OPTIONAL_TARGETS=GetQuadCenters_cli testEdgeMesh ConvertToObj_cli ExpandQuad_cli
ALL_TARGETS=$(DEFAULT_TARGETS) $(OPTIONAL_TARGETS)

CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DBENCHMARK
CPPFLAGS+=-O2 -g
# CPPFLAGS+=-O0 -g
CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP

default_build: $(DEFAULT_TARGETS)
all: $(ALL_TARGETS)

GetQuadDeformations_cli: $(GETQUADDEFS_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

TileQuadFromFile_cli: $(TILEQUAD_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

FindCommonDeformations_cli: $(MAP_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

UnMap_cli: $(UNMAP_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

PatternOptimization_DeformedCell_cli: $(DEFPATOPT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

GetQuadCenters_cli: $(GETQUADCTRS_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

testEdgeMesh: $(testEdgeMesh_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

ExpandQuad_cli: $(EXPQUAD_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

ConvertToObj_cli: $(CONVERT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

%.o: %.cc Makefile
	$(CXX) $(CPPFLAGS) -c $< -o $@

%.o: %.c Makefile
	$(CC) -c $(CFLAGS) $< -o $@

depend:
	@touch Makefile.depend;
	makedepend -Y -f Makefile.depend --  -- -I$(MeshFEM) $(SOURCES) &> /dev/null

clean:
	rm -f $(TARGETS) $(OBJS) *.bak

.PHONY: clean depend

# Read in the dependency file, if it exists
sinclude Makefile.depend

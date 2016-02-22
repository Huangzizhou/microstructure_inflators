include $(MeshFEM)/platform_defs.mk

INCLUDES=-I$(LEVMAR_INC) # bring levmar in first since we have a modified copy that might be overriden by the system's
INCLUDES+=-I$(MeshFEM) -I$(CSGFEM) -isystem $(EIGEN_INC) -I$(SUITESPARSE_INC) -I$(CLIPPER_INC) -I$(TRIANGLE_INC)  \
		  -isystem $(VCGLIB_INC) -isystem $(DLIB_INC)
INCLUDES+=-isystem $(CERES_INC) -isystem $(LIBMATHEVAL_INC)
INCLUDES+=-I$(PYMESH_PATH)/src -I$(PYMESH_PATH)/tools
INCLUDES+=-isystem ../Luigi/wireinflator2D/src

LIBS=$(BOOST_LFLAGS) $(SUITESPARSE_LFLAGS) $(CERES_LFLAGS) $(LIBMATHEVAL_LFLAGS)
LIBS+=$(CLIPPER_LFLAGS) -L$(PYMESH_PATH)/lib -lwires $(TRIANGLE_LFLAGS) 
LIBS+=$(LEVMAR_LFLAGS)


SOURCES=PatternOptimization_cli.cc \
		PatternOptimization_DeformedCell_cli.cc \
		GetQuadDeformations_cli.cc GetQuadCenters_cli.cc \
		TileQuadFromFile_cli.cc ConvertToObj_cli.cc ExpandQuad_cli.cc testEdgeMesh.cc FindCommonDeformations_cli.cc UnMap_cli.cc \
		PatternOptimizationJob.cc RandomJob.cc SpacedJob.cc GradientComponentValidation.cc Inflator_cli.cc ParamRegion_cli.cc
PATOPT_OBJS=PatternOptimization_cli.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o
DEFPATOPT_OBJS=PatternOptimization_DeformedCell_cli.o PatternOptimizationJob.o \
            $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
            $(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o
GETQUADDEFS_OBJS=GetQuadDeformations_cli.o
GETQUADCTRS_OBJS=GetQuadCenters_cli.o
MAP_OBJS=FindCommonDeformations_cli.o UnMap_cli.o
UNMAP_OBJS=UnMap_cli.o
TILEQUAD_OBJS=TileQuadFromFile_cli.o \
            $(MeshFEM)/MeshIO.o \
			$(MeshFEM)/Types.o
testEdgeMesh_OBJS=testEdgeMesh.o \
            $(MeshFEM)/MeshIO.o \
			$(MeshFEM)/Types.o
EXPQUAD_OBJS=ExpandQuad_cli.o \
            $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o \
            $(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o
CONVERT_OBJS=ConvertToObj_cli.o \
            $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o
CREATJOB_OBJS=CreateJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/GlobalBenchmark.o
RANDJOB_OBJS=RandomJob.o  PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/GlobalBenchmark.o
SPACED_OBJS=SpacedJob.o  PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/GlobalBenchmark.o
GRADVALID_OBJS=GradientComponentValidation.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o
INFLATOR_OBJS=Inflator_cli.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o
PREGION_OBJS=ParamRegion_cli.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o
OBJS=$(PATOPT_OBJS) $(DEFPATOPT_OBJS) $(QUADPATOPT_OBJS) $(GIVDEFPATOPT_OBJS) $(CREATJOB_OBJS) $(RANDJOB_OBJS) $(GRADVALID_OBJS) $(SPACED_OBJS) $(INFLATOR_OBJS) $(PREGION_OBJS) $(GETQUADDEFS_OBJS) $(TILEQUAD_OBJS) $(CONVERT_OBJS) $(EXPQUAD_OBJS) $(testEdgeMesh_OBJS) $(GETQUADCTRS_OBJS) $(MAP_OBJS) $(UNMAP_OBJS)

TARGETS=PatternOptimization_cli \
		PatternOptimization_DeformedCell_cli \
		GetQuadDeformations_cli GetQuadCenters_cli \
		TileQuadFromFile_cli testEdgeMesh \
		FindCommonDeformations_cli UnMap_cli \
		CreateJob RandomJob SpacedJob GradientComponentValidation Inflator_cli Inflator2D_cli ParamRegion_cli ConvertToObj_cli ExpandQuad_cli

CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DBENCHMARK
CPPFLAGS+=-O2 -g
# CPPFLAGS+=-O0 -g
CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP

all: $(TARGETS)

PatternOptimization_cli: $(PATOPT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

GetQuadDeformations_cli: $(GETQUADDEFS_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

GetQuadCenters_cli: $(GETQUADCTRS_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

TileQuadFromFile_cli: $(TILEQUAD_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

FindCommonDeformations_cli: $(MAP_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

UnMap_cli: $(UNMAP_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

testEdgeMesh: $(testEdgeMesh_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

ExpandQuad_cli: $(EXPQUAD_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

ConvertToObj_cli: $(CONVERT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

PatternOptimization_DeformedCell_cli: $(DEFPATOPT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

RandomJob: $(RANDJOB_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

SpacedJob: $(SPACED_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

CreateJob: $(CREATJOB_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

GradientComponentValidation: $(GRADVALID_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

Inflator_cli: $(INFLATOR_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

Inflator2D_cli: $(filter-out Inflator_cli.o, $(INFLATOR_OBJS)) Inflator2D_cli.o
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

Inflator2D_cli.o: Inflator_cli.cc Makefile
	$(CXX) -DDIMENSIONS=2 $(CPPFLAGS) -c $< -o $@

ParamRegion_cli: $(PREGION_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

%.o: %.cc Makefile
	$(CXX) $(CPPFLAGS) -c $< -o $@

%.o: %.c Makefile
	$(CC) -c $(CFLAGS) $< -o $@

depend:
	@touch Makefile_m01.depend;
	makedepend -Y -f Makefile_m01.depend --  -- -I$(MeshFEM) $(SOURCES) &> /dev/null

clean:
	rm -f $(TARGETS) $(OBJS) *.bak

.PHONY: clean depend

# Read in the dependency file, if it exists
sinclude Makefile_m01.depend

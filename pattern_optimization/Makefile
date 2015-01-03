include $(MeshFEM)/platform_defs.mk

SOURCES=PatternOptimization_cli.cc PatternOptimizationJob.cc RandomJob.cc SpacedJob.cc GradientComponentValidation.cc
PATOPT_OBJS=PatternOptimization_cli.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o
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
OBJS=$(PATOPT_OBJS) $(CREATJOB_OBJS) $(RANDJOB_OBJS) $(GRADVALID_OBJS) $(SPACED_OBJS) $(INFLATOR_OBJS)

TARGETS=PatternOptimization_cli CreateJob RandomJob SpacedJob GradientComponentValidation Inflator_cli

INCLUDES+=-I$(MeshFEM) -I../Luigi/wireinflator2D/src -isystem $(VCG_INC) -isystem $(DLIB_INC) -I$(PYMESH_PATH)/src -I$(PYMESH_PATH)/tools
LIBS+=-L$(CLIPPER_PATH)/lib -lpolyclipping -L/$(PYMESH_PATH)/lib -lwires
LIBS+=-L$(TRIANGLE_PATH) -ltriangle
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DBENCHMARK
CPPFLAGS+=-O2
# CPPFLAGS+=-O0 -g
CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP

all: $(TARGETS)

PatternOptimization_cli: $(PATOPT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

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

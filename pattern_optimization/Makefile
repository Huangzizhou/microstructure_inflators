MeshFEM ?= $(MESHFEM_PATH)
include $(MeshFEM)/platform_defs.mk

ifeq ($(DEBUG), 1)
	CPPFLAGS+=-O0 -g
else
	CPPFLAGS+=-O2
endif

ifeq ($(DNDEBUG), 1)
	CPPFLAGS+=-DNDEBUG
endif

################################################################################
# Configuration section
################################################################################
# **** Configuration: Use the isosurface inflator for 3D ****
CPPFLAGS+=-DISOSURFACE_INFLATOR -DBENCHMARK
CPPFLAGS+=-Wall -pedantic -std=c++11
# CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP
################################################################################

INCLUDES=-I$(LEVMAR_INC) # bring levmar in first since we have a modified copy that might be overriden by the system's
INCLUDES+=-I$(MeshFEM) -I$(CSGFEM) -I$(EIGEN_INC) -isystem $(SUITESPARSE_INC) -isystem $(BOOST_INC) \
		 -isystem $(CERES_INC) -isystem $(LIBMATHEVAL_INC)
INCLUDES+=-isystem $(CLIPPER_INC) -isystem $(VCGLIB_INC) -isystem $(DLIB_INC)
INCLUDES+=-I$(PYMESH_PATH)/src -I$(PYMESH_PATH)/tools
INCLUDES+=-isystem $(TRIANGLE_INC)

INCLUDES+=-isystem $(ADEPT_INC) -isystem $(VCGLIB_INC) \
		  -isystem $(CGAL_INC) -isystem $(GMP_INC) -isystem $(MPFR_INC)

CPPFLAGS+=$(INCLUDES)

LIBS=$(BOOST_LFLAGS) $(SUITESPARSE_LFLAGS) $(CERES_LFLAGS) $(LIBMATHEVAL_LFLAGS)
LIBS+=$(CLIPPER_LFLAGS) $(TRIANGLE_LFLAGS) $(LEVMAR_LFLAGS)
LIBS+=-L$(PYMESH_PATH)/lib -lwires

ISOINFLATOR_LIBS=-L$(CGAL_LIB) -lCGAL -L$(GMP_LIB) -lgmp -L$(MPFR_LIB) -lmpfr
ISOINFLATOR_LIBS+=-L$(ADEPT_LIB) -ladept
LIBS+=$(ISOINFLATOR_LIBS)

ISOINFLATOR_OBS=../isosurface_inflator/IsosurfaceInflator.o \
				../isosurface_inflator/CGALClippedVolumeMesher.o \
				../isosurface_inflator/VCGSurfaceMesher.o \
				../isosurface_inflator/MeshingOptions.o  \
				../isosurface_inflator/BoxIntersectionMesher.o \
				../isosurface_inflator/MidplaneMesher.o

INFLATOR_OBJS=MakeInflator.o \
				inflators/EqualityConstrainedInflator.o \
				inflators/IsoinflatorWrapper.o \
				inflators/JamesInflatorWrapper.o \
				inflators/LuigiInflatorWrapper.o \
				inflators/BoundaryPerturbationInflator.o

OPTIMIZER_OBJS=optimizers/ceres.o optimizers/dlib.o optimizers/gradient_descent.o \
			   optimizers/nlopt.o

SOURCES=PatternOptimization_cli.cc PatternOptimizationJob.cc RandomJob.cc SpacedJob.cc GradientComponentValidation.cc \
		Inflator_cli.cc ParamRegion_cli.cc
PATOPT_OBJS=PatternOptimization_cli.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o \
			$(INFLATOR_OBJS) $(ISOINFLATOR_OBS) $(OPTIMIZER_OBJS)
CREATJOB_OBJS=CreateJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/GlobalBenchmark.o
RANDJOB_OBJS=RandomJob.o  PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)
SPACED_OBJS=SpacedJob.o  PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)
GRADVALID_OBJS=GradientComponentValidation.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)
INFLATOR_CLI_OBJS=Inflator_cli.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o  $(ISOINFLATOR_OBS)
PREGION_OBJS=ParamRegion_cli.o PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)
OBJS=$(PATOPT_OBJS) $(CREATJOB_OBJS) $(RANDJOB_OBJS) $(GRADVALID_OBJS) $(SPACED_OBJS) $(INFLATOR_CLI_OBJS) $(PREGION_OBJS)

DEFAULT_TARGETS=PatternOptimization_cli # CreateJob GradientComponentValidation Inflator_cli Inflator2D_cli
OPTIONAL_TAGETS=SpacedJob ParamRegion_cli
ALL_TARGETS=$(DEFAULT_TARGETS) $(OPTIONAL_TAGETS)

# DEPENDENCY STUFF
DEPDIR:=.d
COMMA:=,
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS=-MT $@ -MMD -MP -MF $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).Td

default_build: $(DEFAULT_TARGETS)
all: $(ALL_TARGETS)

PatternOptimization_cli: $(PATOPT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@ $(LIBS) # static library ordering

RandomJob: $(RANDJOB_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

SpacedJob: $(SPACED_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

CreateJob: $(CREATJOB_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

GradientComponentValidation: $(GRADVALID_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

Inflator_cli: $(INFLATOR_CLI_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

Inflator2D_cli: $(filter-out Inflator_cli.o, $(INFLATOR_CLI_OBJS)) Inflator2D_cli.o
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

Inflator2D_cli.o: Inflator_cli.cc Makefile
	$(CXX) -DDIMENSIONS=2 $(CPPFLAGS) -c $< -o $@

ParamRegion_cli: $(PREGION_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

# Delete built-in rule, forcing ours to run
%.o: %.cc
%.o: %.cc Makefile
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) -c $< -o $@
	mv -f $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).Td $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).d

%.o: %.c
%.o: %.c Makefile
	$(CC) $(DEPFLAGS) -c $(CFLAGS) $< -o $@
	mv -f $(DEPDIR)/$$(subst .,$(COMMA),$(subst /,__,$*)).Td $(DEPDIR)/$(subst .,$(COMMA),$(subst /,__,$*)).d

clean:
	rm -f $(ALL_TARGETS) $(OBJS) *.bak

.PHONY: clean

# Read dependencies
-include $(DEPDIR)/*.d

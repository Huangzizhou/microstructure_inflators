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
CPPFLAGS+=-DISOSURFACE_INFLATOR -DBENCHMARK -DTOO_LARGE_FOR_METIS
# CPPFLAGS+= -DBENCHMARK
CPPFLAGS+=-Wall -pedantic -Wunused-parameter -Wpedantic -std=c++14
# CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP
################################################################################
INCLUDES=-I$(MeshFEM) -I$(CSGFEM) -I$(EIGEN_INC) -isystem $(SUITESPARSE_INC) -isystem $(BOOST_INC) \
		 -isystem $(CERES_INC) -isystem $(LIBMATHEVAL_INC) -isystem $(CLIPPER_INC) \
		  -I../pattern_optimization -isystem ../Luigi/wireinflator2D/src

INCLUDES+=-I$(PYMESH_PATH)/src -I$(PYMESH_PATH)/tools
INCLUDES+=-isystem $(TRIANGLE_INC)

INCLUDES+=-isystem $(ADEPT_INC) -isystem $(VCGLIB_INC) -isystem $(DLIB_INC) \
		  -isystem $(CGAL_INC) -isystem $(GMP_INC) -isystem $(MPFR_INC)

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

CPPFLAGS+=$(INCLUDES)

M2MSTRESS_OBJS=MacroMicroStress_cli.o  $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o $(MeshFEM)/GlobalBenchmark.o
WCS_OBJS=WCS_cli.o  $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o $(MeshFEM)/GlobalBenchmark.o
GRADVALID_OBJS=GradientComponentValidation.o ../pattern_optimization/PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)
WCSOPT_OBJS=WCSOptimization_cli.o ../pattern_optimization/PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)
WCS_REMESH_OPT_OBJS=WCSRemeshingOptimization_cli.o ../pattern_optimization/PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)

LPHOLE_GRADVALID_OBJS=LpHoleGradientComponentValidation.o ../pattern_optimization/PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)

LPHOLE_INFLATOR_OBJS=LpHoleInflator_cli.o ../pattern_optimization/PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)

OBJS=$(M2MSTRESS_OBJS) $(WCS_OBJS)$(GRADVALID_OBJS) $(WCSOPT_OBJS) $(WCS_REMESH_OPT_OBJS) $(LPHOLE_GRADVALID_OBJS) $(LPHOLE_INFLATOR_OBJS)

TARGETS=MacroMicroStress_cli WCS_cli GradientComponentValidation LpHoleGradientComponentValidation WCSOptimization_cli LpHoleInflator_cli WCSRemeshingOptimization_cli

all: $(TARGETS)

MacroMicroStress_cli: $(M2MSTRESS_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

WCS_cli: $(WCS_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

GradientComponentValidation: $(GRADVALID_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

LpHoleGradientComponentValidation: $(LPHOLE_GRADVALID_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

WCSOptimization_cli: $(WCSOPT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

WCSRemeshingOptimization_cli: $(WCS_REMESH_OPT_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

LpHoleInflator_cli: $(LPHOLE_INFLATOR_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

%.o: %.cc Makefile
	$(CXX) $(CPPFLAGS) -c $< -o $@

%.o: %.c Makefile
	$(CC) -c $(CFLAGS) $< -o $@

depend:
	@touch Makefile.depend;
	makedepend -Y -f Makefile.depend --  -- $(SOURCES) &> /dev/null

clean:
	rm -f $(TARGETS) $(OBJS) *.bak

.PHONY: clean depend

# Read in the dependency file, if it exists
sinclude Makefile.depend

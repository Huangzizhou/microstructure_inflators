MeshFEM ?= $(MESHFEM_PATH)
include $(MeshFEM)/platform_defs.mk

################################################################################
# Configuration section
################################################################################
# **** Configuration: Use the isosurface inflator for 3D ****
CPPFLAGS+=-DISOSURFACE_INFLATOR -DBENCHMARK
CPPFLAGS+=-Wall -pedantic -std=c++11
CPPFLAGS+=-O2
# CPPFLAGS+=-O0 -g
# CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP
################################################################################
INCLUDES+=-I$(MeshFEM) -I$(CSGFEM) -I$(EIGEN_INC) -isystem $(SUITESPARSE_INC) -isystem $(BOOST_INC)
INCLUDES+=-isystem $(CLIPPER_INC) -isystem $(VCGLIB_INC)
INCLUDES+=-I$(PYMESH_PATH)/src -I$(PYMESH_PATH)/tools
INCLUDES+=-isystem $(TRIANGLE_INC)
INCLUDES+=-isystem $(LIBMATHEVAL_INC)
INCLUDES+=-I../../../Luigi/wireinflator2D/src

INCLUDES+=-isystem $(ADEPT_INC) -isystem $(VCGLIB_INC) \
		  -isystem $(CGAL_INC) -isystem $(GMP_INC) -isystem $(MPFR_INC)

CPPFLAGS+=$(INCLUDES)

LIBS=$(BOOST_LFLAGS) $(SUITESPARSE_LFLAGS) $(LIBMATHEVAL_LFLAGS)
LIBS+=$(CLIPPER_LFLAGS) $(TRIANGLE_LFLAGS)
LIBS+=-L$(PYMESH_PATH)/lib -lwires

ISOINFLATOR_LIBS=-L$(CGAL_LIB) -lCGAL -L$(GMP_LIB) -lgmp -L$(MPFR_LIB) -lmpfr
ISOINFLATOR_LIBS+=-L$(ADEPT_LIB) -ladept
LIBS+=$(ISOINFLATOR_LIBS)

ISOINFLATOR_OBS=../../../isosurface_inflator/IsosurfaceInflator.o \
				../../../isosurface_inflator/CGALClippedVolumeMesher.o \
				../../../isosurface_inflator/VCGSurfaceMesher.o \
				../../../isosurface_inflator/MeshingOptions.o  \
				../../../isosurface_inflator/BoxIntersectionMesher.o

SOURCES=PatternOptimization_cli.cc PatternOptimizationJob.cc RandomJob.cc SpacedJob.cc GradientComponentValidation.cc \
		Inflator_cli.cc ParamRegion_cli.cc
SV_VALID_OBJS=ShapeVelocityValidation.o ../../PatternOptimizationJob.o \
			$(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o \
			$(MeshFEM)/EdgeFields.o $(MeshFEM)/GlobalBenchmark.o $(ISOINFLATOR_OBS)
OBJS=$(SV_VALID_OBJS)

TARGETS=ShapeVelocityValidation

all: $(TARGETS)

ShapeVelocityValidation: $(SV_VALID_OBJS)
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

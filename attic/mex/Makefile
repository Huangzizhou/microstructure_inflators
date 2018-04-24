################################################################################
# Objects beginning with "mex" and targets ending with $(MEXEXT) are
# compiled/linked using mex.
#
# The following environment variables should be set:
#	MeshFEM, CSGFEM, MICRO_DIR, TRIANGLE_PATH, CLIPPER_PATH, PYMESH_PATH
# The following platform-specific variables should be configured in
# $MeshFEM/platform_config:
# 	SUITESPARSE_INC, SUITESPARSE_LIBS, BOOST_INC, MATLABDIR, EIGEN_INC,
# 	ARCH
# The following will probably be needed in the future:
# 	VCG_INC, TRIANGLE_LIBS, CLIPPER_LIBS
################################################################################
include $(MeshFEM)/platform_defs.mk

debug     = yes

CXX       ?= g++
CC        ?= gcc
CXXFLAGS += -Wall
LDFLAGS   =

PERHOMO_OBJS=Homogenize.o $(MeshFEM)/Types.o $(MeshFEM)/Materials.o $(MeshFEM)/MeshIO.o
			
INCLUDES+=-I$(MeshFEM) -I$(CSGFEM) -isystem $(SUITESPARSE_INC) -isystem $(BOOST_INC)
INCLUDES+=-I$(MICRO_DIR)/Luigi/wireinflator2D/src -isystem $(EIGEN_INC) # -isystem $(VCG_INC) 

CXXFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES)
CXXFLAGS+=-O2

MEX       ?= $(MATLABDIR)/bin/mex
MEXLIBDIR ?= $(MATLABDIR)/bin/$(ARCH)
MEXEXT    ?= mex$(ARCH)

MEXFLAGS   += -cxx -largeArrayDims -O CXX=$(CXX) -I$(MATLABDIR)/extern/include/ -I$(EIGEN_INC)
MEXLDFLAGS = -L$(MEXLIBDIR)
MEXLDLIBS  += -largeArrayDims $(SUITESPARSE_LIBS) # $(TRIANGLE_LIBS) $(CLIPPER_LIBS)
LDLIBS  += $(SUITESPARSE_LIBS) # $(TRIANGLE_LIBS) $(CLIPPER_LIBS)

ifeq ($(debug), yes)
  MEXFLAGS += -v
endif

OBJS = mex_homogenize.o $(PERHOMO_OBJS)

.SECONDARY:

all: homogenize.$(MEXEXT)

clean:
	-$(RM) *.o *h.gch *~ *.x *.$(MEXEXT)

homogenize.$(MEXEXT): $(OBJS)
	$(MEX) $(MEXLDFLAGS) -output $@ $^ $(MEXLDLIBS)

mex%.o: mex%.cc
	$(MEX) $(MEXFLAGS) -c $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^

MeshFEM ?= $(MESHFEM_PATH)
include $(MeshFEM)/platform_defs.mk

INCLUDES=-I$(MeshFEM) -I$(CSGFEM) -isystem $(EIGEN_INC) -isystem $(BOOST_INC) \
		 -isystem $(ADEPT_INC) -isystem $(VCGLIB_INC) \
		 -isystem $(CGAL_INC) -isystem $(GMP_INC) -isystem $(MPFR_INC) \
		 -isystem $(TRIANGLE_INC)
CPPFLAGS=-O3
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DBENCHMARK

LIBS+=-L$(CGAL_LIB) -lCGAL -L$(GMP_LIB) -lgmp -L$(MPFR_LIB) -lmpfr
LIBS+=$(BOOST_LFLAGS)

all: enumerate enumerate_2D match sort pattern_info test_intersection_check

enumerate: enumerate.o $(MeshFEM)/MeshIO.o intersection_check.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

enumerate_2D: enumerate_2D.o $(MeshFEM)/MeshIO.o intersection_check.o
	g++ -DDIM=2 -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

match: match.o $(MeshFEM)/MeshIO.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

sort: sort.o $(MeshFEM)/MeshIO.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

pattern_info: pattern_info.o $(MeshFEM)/MeshIO.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

test_intersection_check: test_intersection_check.o intersection_check.o $(MeshFEM)/MeshIO.o
	g++ -L/usr/local/lib -L/opt/local/lib $(LIBS) $^ -o $@

enumerate_2D.o: enumerate.cc
	$(CXX) $(CPPFLAGS) -DDIM=2 enumerate.cc -c -o enumerate_2D.o


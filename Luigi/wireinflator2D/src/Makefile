# Creates static libraries for the clipper/triangle dependencies.

INCLUDES=-I$(MeshFEM) -I$(CSGFEM) -isystem $(VCGLIB_PATH) -isystem /opt/local/include/eigen3 -isystem /opt/local/include -isystem$(CLIPPER_PATH)/include/polyclipping
LIBS=-L$(TRIANGLE_PATH) -L$(CLIPPER_PATH)/lib/ -lpolyclipping -ltriangle

CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES)
CPPFLAGS+=-O2
# CPPFLAGS+=-g -O0

all: inflator

inflator: main.cpp
	$(CXX) $(CPPFLAGS) main.cpp $(MeshFEM)/MeshIO.cc $(MeshFEM)/EdgeFields.cc $(LIBS) -o inflator

clean:
	rm inflator

.PHONY: clean depend

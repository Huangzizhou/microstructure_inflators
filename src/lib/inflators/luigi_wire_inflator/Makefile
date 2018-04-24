# Creates static libraries for the clipper/triangle dependencies.
all: libtriangle.a libclipper.a

CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES)
CPPFLAGS+=-O2

build_dir:
	mkdir -p build

lib_dir:
	mkdir -p lib

libtriangle.a: src/triangle.c src/triangle.h build_dir lib_dir
	$(CC) $(CFLAGS) -c src/triangle.c -o build/triangle.o
	ar rvs lib/libtriangle.a build/triangle.o
libclipper.a: src/clipper.cpp src/clipper.hpp build_dir
	$(CXX) $(CPPFLAGS) -c src/clipper.cpp -o build/clipper.o
	ar rvs lib/libclipper.a build/clipper.o

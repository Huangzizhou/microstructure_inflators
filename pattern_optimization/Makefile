include $(MeshFEM)/platform_defs.mk

PATOPT_OBJS=PatternOptimization_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o \
			$(MeshFEM)/Materials.o $(MeshFEM)/MSHFieldParser.o
OBJS=$(PATOPT_OBJS)

TARGETS=PatternOptimization_cli

INCLUDES+=-I$(MeshFEM) -I../Luigi/wireinflator2D/src -isystem $(HOME)/Research/vcglib/
LIBS+=../Luigi/wireinflator2D/lib/libclipper.a
LIBS+=../Luigi/wireinflator2D/lib/libtriangle.a
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES)
CPPFLAGS+=-O2
# CPPFLAGS+=-O0 -g
CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP

all: $(TARGETS)

PatternOptimization_cli: $(PATOPT_OBJS)
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

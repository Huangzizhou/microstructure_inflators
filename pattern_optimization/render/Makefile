include $(MeshFEM)/platform_defs.mk

RENDER_OBJS=Render_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/EdgeFields.o \
			$(CSGFEM)/draw.o
OBJS=$(RENDEROPT_OBJS)
TARGETS=Render_cli

INCLUDES+=-I$(MeshFEM) -I/opt/local/include/freetype2
LIBS+=-losmesa -lpng
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DUSE_MESA
CPPFLAGS+=-O2
# CPPFLAGS+=-O0 -g
CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP

all: $(TARGETS)

Render_cli: $(RENDER_OBJS)
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

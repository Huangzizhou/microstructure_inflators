include $(MeshFEM)/platform_defs.mk

RENDER_OBJS=Render_cli.o $(MeshFEM)/MeshIO.o $(MeshFEM)/Types.o \
			$(MeshFEM)/MSHFieldParser.o $(MeshFEM)/EdgeFields.o \
			$(CSGFEM)/draw.o
IVIEW_OBJS=inflate_viewer.o $(MeshFEM)/MeshIO.o 
IVIEW_LIBS=-framework GLUT -framework OpenGL -L$(CLIPPER_PATH)/lib -lpolyclipping -L/opt/local/lib -lAntTweakBar -lmatheval -ltriangle
OBJS=$(RENDEROPT_OBJS) $(IVIEW_OBJS)
TARGETS=Render_cli inflate_viewer

INCLUDES+=-I$(MeshFEM) -I/opt/local/include/freetype2 -I$(CLIPPER_PATH)/include -I$(EIGEN_INC) -I.. -I../../Luigi/wireinflator2D/src -isystem $(VCG_INC) -I$(PYMESH_PATH)/src -I$(PYMESH_PATH)/tools
LIBS+=$(RENDER_LIBS)
CPPFLAGS+=-Wall -pedantic -std=c++11 $(INCLUDES) -DUSE_MESA
CPPFLAGS+=-O2
# CPPFLAGS+=-O0 -g
CPPFLAGS+=-DHAVE_NAMESPACES -DHAVE_STD # Garbage for OptPP

all: $(TARGETS)

Render_cli: $(RENDER_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(LIBS) -o $@

inflate_viewer: $(IVIEW_OBJS)
	$(CXX) $(CPPFLAGS) $^ $(IVIEW_LIBS) -o $@

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

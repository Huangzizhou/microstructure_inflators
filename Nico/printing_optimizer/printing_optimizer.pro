macx: QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9
QMAKE_MAC_SDK = macosx10.9
CONFIG += c++11

VCGLIBDIR = ../../../vcg/vcglib/
GLEWDIR   = ../../../code/lib/glew/
ANTDIR    = ../../../code/lib/AntTweakBar1.16

HEADERS       = glwidget.h \
                edge_mesh_type.h\
                tri_mesh_type.h

SOURCES       = glwidget.cpp \
                main.cpp \

QT           += opengl

# Compile glew
DEFINES += GLEW_STATIC
DEFINES += USE_OPENGL
INCLUDEPATH += $$GLEWDIR/include
SOURCES += $$GLEWDIR/src/glew.c

INCLUDEPATH += $$VCGLIBDIR
INCLUDEPATH += $$GLEWDIR/include
INCLUDEPATH += $$ANTDIR/include

SOURCES += $$VCGLIBDIR/wrap/ply/plylib.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackball.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackmode.cpp
SOURCES += $$VCGLIBDIR/wrap/qt/anttweakbarMapperNew.cpp

# Awful problem with windows..
win32{
  DEFINES += NOMINMAX
  LIBS +=$$ANTDIR/lib/AntTweakBar.lib
}

mac {
    # Mac specific config required to avoid to make application bundles
    CONFIG -= app_bundle

    ANTTWEAKLIB = libAntTweakBar.dylib

   LIBS += $$ANTDIR/lib/$$ANTTWEAKLIB
   QMAKE_POST_LINK += "cp -P $$ANTDIR/lib/$$ANTTWEAKLIB ./libAntTweakBar.dylib ;"
   QMAKE_POST_LINK += "install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib ./$$TARGET ;"
}

unix:!macx{
  LIBS +=$$ANTDIR/lib/libAntTweakBar.so -lGLU
}

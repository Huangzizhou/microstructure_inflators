macx: QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9
QMAKE_MAC_SDK = macosx10.9
CONFIG += c++11

VCGLIBDIR = /Users/qingnanzhou/Research/sources/vcglib

INCLUDEPATH += $$VCGLIBDIR

HEADERS       = edge_mesh_type.h\
                tri_mesh_type.h

SOURCES       = main_command.cpp 


mac {
    # Mac specific config required to avoid to make application bundles
    CONFIG -= app_bundle
}

unix:!macx{
  LIBS +=$$ANTDIR/lib/libAntTweakBar.so -lGLU
}

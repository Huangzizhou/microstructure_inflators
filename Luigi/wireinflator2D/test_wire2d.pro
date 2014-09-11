VCGLIBDIR   = ../../../vcglib

DESTDIR     = ./bin
BUILD_DIR   = ./bin
UI_DIR      = $$BUILD_DIR
MOC_DIR     = $$BUILD_DIR
OBJECTS_DIR = $$BUILD_DIR
RCC_DIR     = $$BUILD_DIR

QT       += core
QT       -= gui

TARGET = test_wire2d

# Enable C++11
CONFIG += c++11

CONFIG -= app_bundle

INCLUDEPATH += $$VCGLIBDIR

HEADERS += \
    src/clipper.hpp \
    src/clipperHelper.h \
    src/EdgeMeshPattern.h \
    src/EdgeMeshType.h \
    src/InflatorParameters.h \
    src/OutMesh.h \
    src/Pattern2D.h \
    src/tessellator2d.h \
    src/triangle.h \
    src/TriMeshType.h \
    src/WireInflator2D.h \
    src/EdgeMeshUtils.h \
    src/WireMesh2D.h

SOURCES += \
    src/clipper.cpp \
    src/triangle.c \
    src/main.cpp

OTHER_FILES += \
	meshes/octa_cell.obj

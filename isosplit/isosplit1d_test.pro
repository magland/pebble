QT       += core
QT       -= gui

TARGET = isosplit1d_test
DESTDIR = bin
OBJECTS_DIR = build
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

SOURCES += isosplit1d_test.cpp

INCLUDEPATH += ../mdaio
DEPENDPATH += ../mdaio
VPATH += ../mdaio
HEADERS += mda.h mdaio.h usagetracking.h
SOURCES += mda.cpp mdaio.cpp usagetracking.cpp

HEADERS += jisotonic.h isosplit1d.h isosplit.h
SOURCES += jisotonic.cpp isosplit1d.cpp isosplit.cpp

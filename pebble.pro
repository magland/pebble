#-------------------------------------------------
#
# Project created by QtCreator 2015-08-28T17:20:25
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = pebble
DESTDIR = bin
OBJECTS_DIR = build
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    pebble_helpers.cpp

HEADERS += parse_command_line_params.h \
    pebble_helpers.h
SOURCES += parse_command_line_params.cpp

HEADERS += do_pca.h
SOURCES += do_pca.cpp

INCLUDEPATH += mdaio
DEPENDPATH += mdaio
VPATH += mdaio
HEADERS += mda.h mdaio.h usagetracking.h
SOURCES += mda.cpp mdaio.cpp usagetracking.cpp

INCLUDEPATH += procedures
DEPENDPATH += procedures
VPATH += procedures
HEADERS += find_critical_times.h
SOURCES += find_critical_times.cpp

INCLUDEPATH += isosplit
DEPENDPATH += isosplit
VPATH += isosplit
HEADERS += jisotonic.h isosplit1d.h isosplit.h
SOURCES += jisotonic.cpp isosplit1d.cpp isosplit.cpp

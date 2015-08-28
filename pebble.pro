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


SOURCES += main.cpp

INCLUDEPATH += mdaio
DEPENDPATH += mdaio
VPATH += mdaio
HEADERS += mda.h mdaio.h usagetracking.h
SOURCES += mda.cpp mdaio.cpp usagetracking.cpp

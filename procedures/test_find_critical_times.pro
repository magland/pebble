QT       += core
QT       -= gui

TARGET = test_find_critical_times
DESTDIR = .
OBJECTS_DIR = build
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

SOURCES += test_find_critical_times.cpp

HEADERS += find_critical_times.h
SOURCES += find_critical_times.cpp

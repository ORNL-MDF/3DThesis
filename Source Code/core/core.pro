#-------------------------------------------------
#
# Project created by QtCreator 2018-10-29T10:43:26
#
#-------------------------------------------------

QT       -= core gui

TARGET = HT_Model_Lib
TEMPLATE = lib
CONFIG += staticlib


INCLUDEPATH += include

QMAKE_CXXFLAGS += -openmp -pthread
LIBS += -openmp

unix {
    target.path = /usr/lib
    INSTALLS += target
}

SOURCES += \
    src/util.cpp \
    src/ht_point.cpp \
    src/driver.cpp

HEADERS += \
    include/data_structs.h \
    include/util.h \
    include/ht_point.h \
    include/driver.h

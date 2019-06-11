#-------------------------------------------------
#
# Project created by QtCreator 2018-11-04T15:28:54
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ../SAHTMpublicRelease
#TARGET = SAHTMpublicRelease
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp \
        modelwrapper.cpp

HEADERS += \
        mainwindow.h \
        modelwrapper.h \
    progresscallback.h

FORMS += \
        mainwindow.ui

QMAKE_CXXFLAGS += -openmp -pthread
LIBS += -openmp

INCLUDEPATH += ../core/include/

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../core/release/ -lHT_Model_Lib
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../core/debug/ -lHT_Model_Lib
else:unix: LIBS += -L$$OUT_PWD/../SAHTM_lib/ -lHT_Model_Lib

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../core/release/libHT_Model_Lib.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../core/debug/libHT_Model_Lib.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../core/release/HT_Model_Lib.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../core/debug/HT_Model_Lib.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../core/libHT_Model_Lib.a



# Default rules for deployment.



qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target




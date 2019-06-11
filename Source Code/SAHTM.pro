#-------------------------------------------------
#
# Project created by QtCreator 2018-11-5T10:51:20
#
#-------------------------------------------------

TEMPLATE = subdirs
CONFIG += ordered
SUBDIRS = core \
          GUI

GUI.depends = core

cmake_minimum_required(VERSION 3.6.0)

project(SAHTM_GUI)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
find_package(Qt5 COMPONENTS Core Widgets REQUIRED)

set(SAHTM_GUI_SRCS
    mainwindow.ui
    mainwindow.cpp
    modelwrapper.cpp
	main.cpp
	progresscallback.h
)

add_executable (SAHTM_GUI ${SAHTM_GUI_SRCS})

target_link_libraries(SAHTM_GUI LINK_PUBLIC SAHTM_Lib Qt5::Core Qt5::Widgets)

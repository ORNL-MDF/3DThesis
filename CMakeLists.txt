cmake_minimum_required(VERSION 3.9)

project(Thesis CXX)

include(GNUInstallDirs)

find_package(OpenMP)

add_subdirectory(src)

cmake_minimum_required(VERSION 3.9)

project(Thesis CXX)

include(GNUInstallDirs)

find_package(OpenMP)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_subdirectory(src)
add_subdirectory(bin)

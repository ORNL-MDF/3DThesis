cmake_minimum_required(VERSION 3.9)

project(Thesis CXX)
set(PROJECT_VERSION "3.1.0-dev")

include(GNUInstallDirs)

find_package(OpenMP)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_subdirectory(src)
add_subdirectory(bin)

cmake_minimum_required(VERSION 3.9)

project(Thesis CXX)
set(PROJECT_VERSION "3.1.0-dev")

include(GNUInstallDirs)

find_package(OpenMP)

if(Thesis_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

add_subdirectory(src)
add_subdirectory(bin)

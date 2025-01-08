cmake_minimum_required(VERSION 3.9)

project(Thesis CXX)
set(PROJECT_VERSION "3.1.0-dev")

include(GNUInstallDirs)

find_package(OpenMP)

find_package(MPI QUIET COMPONENTS CXX)
option(Thesis_REQUIRE_MPI "Require Thesis to build with MPI support" MPI_FOUND)
if(Thesis_REQUIRE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
set(Thesis_ENABLE_MPI ${MPI_FOUND})

add_subdirectory(src)
add_subdirectory(bin)

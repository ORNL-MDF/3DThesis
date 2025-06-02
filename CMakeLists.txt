cmake_minimum_required(VERSION 3.9)

project(Thesis CXX)
set(PROJECT_VERSION "4.1.0-dev")

include(GNUInstallDirs)

find_package(OpenMP)

find_package(MPI QUIET COMPONENTS CXX)
option(Thesis_REQUIRE_MPI "Require Thesis to build with MPI support" MPI_FOUND)
if(Thesis_REQUIRE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
set(Thesis_ENABLE_MPI ${MPI_FOUND})

find_package(Git)
if(GIT_FOUND AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log --pretty=format:%H -n 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE Thesis_COMMIT_HASH)
else()
  set(Thesis_COMMIT_HASH "Not a git repository")
endif()
message(STATUS "Thesis version: ${PROJECT_VERSION}")
message(STATUS "Thesis commit:  ${Thesis_COMMIT_HASH}")

add_subdirectory(src)
add_subdirectory(bin)

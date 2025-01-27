cmake_minimum_required(VERSION 3.16)
project(Thesis LANGUAGES CXX)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MPI if enabled
if(THESIS_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

# Find Kokkos
find_package(Kokkos REQUIRED)

# Find Stork
find_package(Stork REQUIRED)

# Gather all .cpp files in the include/ directory
file(GLOB_RECURSE SRC_FILES include/*.cpp)

# Gather all .h/.hpp files in the include/ directory
file(GLOB_RECURSE INCLUDE_FILES include/*.hpp)

# Generate ThesisConfig.h
configure_file(
    include/ThesisConfig.h.cmake.in
    ${CMAKE_BINARY_DIR}/include/ThesisConfig.h
)

# Core library
add_library(Thesis ${SRC_FILES})

# Link the dependencies to the Thesis library
target_link_libraries(Thesis
    PUBLIC
    Stork::Stork
    Kokkos::kokkos
)

# Link OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  target_link_libraries(Thesis PUBLIC OpenMP::OpenMP_CXX)
endif()

# Conditionally link MPI
if(THESIS_ENABLE_MPI)
  target_link_libraries(Thesis PUBLIC MPI::MPI_CXX)
endif()

# Include directories for Thesis
target_include_directories(Thesis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>                 # Generated ThesisConfig.h
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  # Public headers
    $<INSTALL_INTERFACE:include>                          # Installation headers
)

# Set the output directory for the executables during build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps)

# Add executables for each file in apps/
file(GLOB APP_SOURCES "apps/*.cpp")  # Find all .cpp files in apps/
foreach(app_source ${APP_SOURCES})
    get_filename_component(app_name ${app_source} NAME_WE)  # Get the file name without extension
    add_executable(${app_name} ${app_source})              # Create an executable for each app
    target_link_libraries(${app_name} PRIVATE Thesis)      # Link each app with Thesis
endforeach()

# Install the Thesis library
install(TARGETS Thesis
    EXPORT ThesisConfig
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install public headers
install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/include/ThesisConfig.h
    DESTINATION include
)

# Package configuration for find_package() support
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ThesisConfigVersion.cmake"
    VERSION 1.0
    COMPATIBILITY AnyNewerVersion
)

# Create the ThesisConfig.cmake file
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/ThesisConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ThesisConfig.cmake"
    INSTALL_DESTINATION lib/cmake/Thesis
)

# Install the export configuration for package support
install(EXPORT ThesisConfig
    DESTINATION lib/cmake/Thesis
    NAMESPACE Thesis::
)

# Install the package configuration files
install(EXPORT ThesisConfig
    FILE ThesisTargets.cmake
    NAMESPACE Thesis::
    DESTINATION lib/cmake/Thesis
)

# Install version and config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ThesisConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ThesisConfigVersion.cmake"
    DESTINATION lib/cmake/Thesis
)

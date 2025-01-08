#!/bin/bash

# Clean previous build
rm -rf build
clear

# Set environment variables
export NVCC_WRAPPER_DEFAULT_COMPILER=mpic++

# Set environment variable for CodeA installation path
export STORK_DIR=$HOME/GitCode/stork/build/install

# Set environment variable for building with MPI
export THESIS_ENABLE_MPI=true

# Create build directory
mkdir -p build
pushd build

# CMake configuration
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_INSTALL_PREFIX=install \
  -D CMAKE_CXX_FLAGS="-fopenmp -O3 -ffast-math -march=znver3 -mtune=znver3" \
  -D CMAKE_PREFIX_PATH="${KOKKOS_INSTALL};${MPI_DIR};${STORK_DIR}" \
  -D CMAKE_CUDA_ARCHITECTURES="86" \
  -D THESIS_ENABLE_MPI=ON \
  ..

# Build the project
make -j install  # Use all available CPU cores for faster compilation

# Return to original directory
popd

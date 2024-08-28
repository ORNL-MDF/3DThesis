
## Dependencies
3DThesis has the following dependencies:

|Dependency | Version  | Required | Details|
|---------- | -------  |--------  |------- |
|CMake      | 3.9+     | No      | Build system
|OpenMP     | -        | Yes      | Parallel threading

## Build and install

### CMake

The following script will configure, build, and install 3DThesis. The only flag passed to CMake creates a local install:
```
cd 3DThesis
mkdir build
cd build
cmake \
    -D CMAKE_INSTALL_PREFIX=install \
    ..
make -j install
```

The following example adds more optional CMake flags to further customize the build:
```
cd 3DThesis
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE="Release" \
    -D CMAKE_CXX_COMPILER=g++ \
    -D CMAKE_CXX_FLAGS="-Wall -Wextra -pedantic" \
    -D CMAKE_INSTALL_PREFIX=install \
    -D CMAKE_PREFIX_PATH="$OMP_INSTALL" \
    ..
make -j install
```

### make

An example makefile is also supported. To perform an in-source build, simply type `make`
from the 3DThesis directory. Compiler settings, build flags, and build location can all
be changed within that file.

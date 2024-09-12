rm -r build
clear
mkdir build
pushd build
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_CXX_FLAGS="-O3 -ffast-math" \
  -D CMAKE_PREFIX_PATH="$KOKKOS_DIR/build/install" \
  -D CMAKE_INSTALL_PREFIX=install \
  .. ;
make -j install
popd

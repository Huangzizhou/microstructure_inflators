#! /bin/bash

# Prepare external libraries (Ceres ...)
pushd 3rdparty
mkdir build
cd build
cmake ..
make -j4
popd

# Prepare build folders
mkdir build
for CONFIG in "Debug" "Release"; do
	mkdir build/${CONFIG}
	pushd build/${CONFIG}
	cmake -DCMAKE_BUILD_TYPE=${CONFIG} ../..
	popd
done

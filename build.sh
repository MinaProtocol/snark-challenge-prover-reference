#!/bin/bash
mkdir build
pushd build
  cmake -DMULTICORE=ON -DUSE_PT_COMPRESSION=OFF $EXTRA_CMAKE_ARGS_FOR_CI ..
  make -j12 main opencl_main generate_parameters
popd
mv build/libsnark/main .
mv build/libsnark/opencl_main .
mv build/libsnark/generate_parameters .

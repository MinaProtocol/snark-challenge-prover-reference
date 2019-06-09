#!/bin/bash
mkdir build
pushd build
  cmake -DMULTICORE=ON -DUSE_PT_COMPRESSION=OFF .. 
  make -j12 main generate_inputs generate_parameters
popd
mv build/libsnark/main .
mv build/libsnark/generate_inputs .
mv build/libsnark/generate_parameters .

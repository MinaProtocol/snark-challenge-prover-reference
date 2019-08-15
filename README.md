### Groth16 prover
This directory contains a reference CPU and GPU implementation of  the
Groth16 prover
using [libsnark](README-libsnark.md).

There are two provers implemented in this repository. The first is `libsnark/main.cpp`, which has _the_ CPU libsnark reference prover that we compare benchmarks against. The second is `libsnark/opencl_main.cpp`.

The GPU prover implements a custom DIT/DIF FFT and Pippenger Multiexponentiation algorithms in OpenCL. The Kernel is currently unoptimized and this is a WIP.

NOTE: Due to an AMD compiler bug we are unable to use this kernel on AMD hardware. We will continue to explore our options with Kronos Group to open this implementation up to both hardware vendors. 

#### Dependencies

The code should compile and run on Ubuntu 18.04 with the following dependencies installed:

``` bash
sudo apt-get install -y build-essential \
    cmake \
    git \
    libomp-dev \
    libgmp3-dev \
    libprocps-dev \
    python-markdown \
    libboost-all-dev \
    libssl-dev \
    pkg-config \
    nvidia-cuda-toolkit
```


Building on MacOS is not recommended as CUDA support is harder to use. (Apple mostly ships with AMD.)


#### Build
``` bash
./build.sh
```

#### Generate parameters and inputs
``` bash
./generate_parameters
```

When iterating on your implementation for correctness (and not performance)
smaller constraint systems are fine. In that case, `./generate_parameters fast`
will give you smaller parameters that don't take as long to generate or
prove with.

### Run
``` bash
./main MNT4753 compute MNT4753-parameters MNT4753-input MNT4753-output
./main MNT6753 compute MNT6753-parameters MNT6753-input MNT6753-output
./opencl_main MNT4753 compute MNT4753-parameters MNT4753-input MNT4753-output
```

### Check results
``` bash
sha256sum MNT4753-output MNT6753-output MNT4753-output_cuda MNT6753-output_cuda
```

‡ Due to some bug with `nvcc` in version `10.1`, including any of the `mnt4`/`mnt6` `libff` headers (a library libsnark uses for finite field arithmetic) leads to a compilation failure. The wrapper library uses something like the [pImpl idiom](https://en.cppreference.com/w/cpp/language/pimpl) to hide the need for those headers from the caller.


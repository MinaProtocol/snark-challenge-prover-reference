### Groth16 prover
This directory contains a reference CPU implementation of  the
Groth16 prover
using [libsnark](README-libsnark.md).

The build system is already set up to integrate CUDA. `cuda-fixnum/main.cu`
will be built and linked against `libsnark/main.cu`. Currently only one
wrapper function is exposed, `do_fixnum_example`. Define the CUDA functions you
need in `cuda-fixnum/main.cu`, expose wrapper functions for them, and call them
as needed from the prover code.

#### Dependancies

The code should compile and run on Ubuntu 18.04 with the following dependencies installed.
You will also need the CUDA toolkit.

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
    pkg-config
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

### Run
``` bash
./main MNT4753 compute CPU MNT4753-parameters MNT4753-input MNT4753-output
./main MNT6753 compute CPU MNT6753-parameters MNT6753-input MNT6753-output
./main MNT4753 compute GPU ./cuda-fixnum/inputs gpu-output
```

### Check results
``` bash
sha256sum MNT4753-output MNT6753-output
```

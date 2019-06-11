### Groth16 prover
This directory contains a reference CPU implementation of  the
Groth16 prover
using [libsnark](README-libsnark.md).

#### Dependancies

The code should compile and run on Ubuntu 18.04 with the following dependancies installed:

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


Building on MacOS is not reccomended as CUDA support is harder to use. (Apple mostly ships with AMD.)


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
./main MNT4753 compute MNT4753-parameters MNT4753-input MNT4753-output
./main MNT6753 compute MNT6753-parameters MNT6753-input MNT6753-output
```

### Check results
``` bash
sha256sum MNT4753-output MNT6753-output
```

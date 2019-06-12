#!/usr/bin/env bash

set -e

sleep 5

#=================================#
#  Download Repositories for Dev  #
#=================================#

# git clone https://github.com/CodaProtocol/snark-challenge-prover-reference.git
git clone https://github.com/CodaProtocol/cuda-fixnum.git

#=================================#
#  Download Dependencies for Dev  #
#=================================#
sudo apt-get update -y 
sudo apt-get install -y build-essential \
    cmake \
    git \
    libgmp3-dev \
    libprocps-dev \
    python-markdown \
    libboost-all-dev \
    libssl-dev 

#=================================#
#     Set GCC & G++ Versions      #
#=================================#
sudo add-apt-repository ppa:ubuntu-toolchain-r/ppa
sudo apt-get update
sudo apt-get install -y gcc-7 g++-7
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 99 --slave /usr/bin/g++ g++ /usr/bin/g++-7

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
    dkms freeglut3 freeglut3-dev libxi-dev libxmu-dev \
    cmake \
    git \
    libgmp3-dev \
    libprocps-dev \
    python-markdown \
    libboost-all-dev \
    libssl-dev 

sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub

wget http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-repo-ubuntu1804_10.1.168-1_amd64.deb
sudo dpkg -i ./cuda-repo-ubuntu1804_10.1.168-1_amd64.deb
sudo apt-get update
sudo apt-get install cuda-10-1 -y
nvidia-smi -pm 1
# Verify Install
nvidia-smi

# Add Cuda Tools to PATH
touch ~/.bash_profile
echo "PATH=/usr/local/cuda-10.1/bin:$PATH" > ~/.bash_profile

# Cleanup
rm cuda-repo-ubuntu1804_10.1.168-1_amd64.deb

#=================================#
#     Set GCC & G++ Versions      #
#=================================#
sudo add-apt-repository ppa:ubuntu-toolchain-r/ppa
sudo apt-get update
sudo apt-get install -y gcc-7 g++-7
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 99 --slave /usr/bin/g++ g++ /usr/bin/g++-7


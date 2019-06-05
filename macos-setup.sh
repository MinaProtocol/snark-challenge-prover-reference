#!/bin/bash
set -x #echo on
set -eu

PACKAGES="gpatch opam cmake gmp pkg-config openssl libffi libsodium boost zlib libomp"

# removing already installed packages from the list
for p in $(env HOMEBREW_NO_AUTO_UPDATE=1 brew list); do
  PACKAGES=${PACKAGES//$p/}
done;

# only run if there's work to do
if [[ $PACKAGES = *[![:space:]]* ]]; then
  yes | env HOMEBREW_NO_AUTO_UPDATE=1 brew install $PACKAGES
else
  echo 'All brew packages have already been installed.'
fi

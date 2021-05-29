#! /bin/bash

source script/discover_os
discover_os

# check if the directory $1/BLITZ exist

if [ -d "$1/BLITZ" ]; then
  echo "BLITZ is already installed"
else
  ## Remove old download
  rm blitz-1.0.2.tar.gz
  rm -rf blitz-1.0.2.tar.gz
  wget https://github.com/blitzpp/blitz/archive/refs/tags/1.0.2.tar.gz -O blitz-1.0.2.tar.gz
  tar -xf blitz-1.0.2.tar.gz
  cd blitz-1.0.2

  BUILDDIR=build
  mkdir -p $BUILDDIR
  cd $BUILDDIR
  echo "cmake ../. -DCMAKE_INSTALL_PREFIX=$1/BLITZ"
  cmake ../. -DCMAKE_INSTALL_PREFIX=$1/BLITZ
  make -j $2
  make install

  # Mark the installation
  echo 1 > $1/BLITZ/version
fi

## Algoim installation




if [ -d "$1/ALGOIM" ]; then
  echo "ALGOIM is already installed"
else

  ## Remove old download
  rm algoim.tar.gz
  rm -rf algoim.tar.gz
  wget http://ppmcore.mpi-cbg.de/upload/algoim.tar.gz
  tar -xf algoim.tar.gz
  mv algoim $1/ALGOIM
  mv $1/ALGOIM/src $1/ALGOIM/include
  # Mark the installation
  echo 1 > $1/ALGOIM/version
fi

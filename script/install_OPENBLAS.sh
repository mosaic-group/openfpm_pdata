#! /bin/bash

# check if the directory $1/OPENBLAS exist

if [ -d "$1/OPENBLAS" ]; then
  echo "OPENBLAS already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/OpenBLAS-0.2.15.tar.gz
rm -rf OpenBLAS-0.2.15
tar -xf OpenBLAS-0.2.15.tar.gz
cd OpenBLAS-0.2.15

# configuration

make CC=gcc CXX=g++
mkdir $1/OPENBLAS
make install PREFIX=$1/OPENBLAS
rm -rf OpenBLAS-0.2.15


#! /bin/bash

# check if the directory $1/HDF5 exist

if [ -d "$1/HDF5" ]; then
  echo "HDF5 already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/hdf5-1.8.16.tar.bz2
tar -xf hdf5-1.8.16.tar.bz2
cd hdf5-1.8.16
CC=mpicc ./configure --enable-parallel --prefix=$1/HDF5
make -j 4
mkdir $1/HDF5
make install

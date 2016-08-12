#! /bin/bash

# check if the directory $1/HDF5 exist

if [ -d "$1/LIBHILBERT" ]; then
  echo "LIBHILBERT already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/libhilbert-master.tar.gz
rm -rf libhilbert-master
tar -xf libhilbert-master.tar.gz
cd libhilbert-master
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/LIBHILBERT ..
make all
make install


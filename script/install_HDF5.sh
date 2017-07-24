#! /bin/bash

source script/discover_os

discover_os

# check if the directory $1/HDF5 exist

if [ -d "$1/HDF5" ]; then
  echo "HDF5 is already installed"
  exit 0
fi

if [ ! -d "$1/ZLIB"  -a x"$platform" != x"cygwin" ]; then
  rm zlib1211.tar.gz
  rm -rf zlib-1.2.11
  wget https://zlib.net/zlib-1.2.11.tar.gz
  if [ $? -ne 0 ]; then
    echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
    exit 1
  fi
  tar -xf zlib1211.tar.gz
  cd zlib-1.2.11

  CC=mpicc ./configure --prefix=$1/ZLIB
  make -j $2

  if [ $? -eq 0 ]; then
    make check install
  else
    echo -e "\033[91;5;1m ZLIB Installation FAILED \033[0m"
    exit 1
  fi

else
  echo "ZLIB is already installed"
fi

wget http://ppmcore.mpi-cbg.de/upload/hdf5-1.8.19.tar.gz
tar -xf hdf5-1.8.19.tar.gz
cd hdf5-1.8.19

if [ x"$plaform" != "cygwin" ]; then
        CC=mpicc ./configure --with-zlib=$1/ZLIB --enable-parallel --prefix=$1/HDF5
else
        CC=mpicc ./configure --enable-parallel --prefix=$1/HDF5
fi
make -j $2
mkdir $1/HDF5
make install
echo 1 > $1/HDF5/version

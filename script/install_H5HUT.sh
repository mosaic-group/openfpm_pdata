#! /bin/bash

# check if the directory $1/H5HUT exist

if [ -d "$1/H5HUT" ]; then
  echo "H5HUT already installed"
  exit 0
fi

rm H5hut-1.99.13.tar.gz
rm -rf H5hut-1.99.13
wget http://ppmcore.mpi-cbg.de/upload/H5hut-1.99.13.tar.gz
tar -xf H5hut-1.99.13.tar.gz
cd H5hut-1.99.13
sh ./autogen.sh
./configure --with-hdf5=$1/HDF5 --enable-parallel --prefix=$1/H5HUT
make -j 4
if [ $? -eq 0 ]; then
    echo "H5HUT error installing"
    exit 0
fi
mkdir $1/H5HUT
make install

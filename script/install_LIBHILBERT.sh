#! /bin/bash

rm -rf libhilbert-master
wget http://ppmcore.mpi-cbg.de/upload/libhilbert-master.tar.gz -O boost_1_84_0.tar.gz
tar -xf libhilbert-master.tar.gz
cd libhilbert-master
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/LIBHILBERT ..
make all -j $2
make install


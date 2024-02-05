#! /bin/bash

rm -rf Vc-1.4.1
wget http://ppmcore.mpi-cbg.de/upload/Vc-1.4.1.tar.gz -O Vc-1.4.1.tar.gz
tar -xf Vc-1.4.1.tar.gz
cd Vc-1.4.1
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/VCDEVEL -DCMAKE_C_COMPILER=$3 -DCMAKE_CXX_COMPILER=$4 ..
make -j $2
make install


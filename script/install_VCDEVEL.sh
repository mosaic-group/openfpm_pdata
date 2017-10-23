#! /bin/bash

# check if the directory $1/VCDEVEL exist

if [ -d "$1/VCDEVEL" ]; then
  echo "VCDEVEL already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/Vc-1.3.2.tar.gz
#rm -rf Vc
tar -xf Vc-1.3.2.tar.gz
cd Vc-1.3.2
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/VCDEVEL ..
make
make install


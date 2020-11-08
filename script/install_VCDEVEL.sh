#! /bin/bash

# check if the directory $1/VCDEVEL exist

if [ -d "$1/VCDEVEL" -a -f "$1/VCDEVEL/include/Vc/Vc" ]; then
  echo "VCDEVEL already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/Vc-1.4.1.tar.gz
#rm -rf Vc
tar -xf Vc-1.4.1.tar.gz
cd Vc-1.4.1
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/VCDEVEL -DCMAKE_C_COMPILER=$3 -DCMAKE_CXX_COMPILER=$4 ..
make
make install


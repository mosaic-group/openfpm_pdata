#! /bin/bash

# check if the directory $1/METIS exist

if [ -d "$1/METIS" ]; then
  echo "METIS already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/metis-5.1.0.tar.gz
tar -xf metis-5.1.0.tar.gz
cd metis-5.1.0
cputype=$(uname -m | sed "s/\\ /_/g")
systype=$(uname -s)
BUILDDIR=build/$systype-$cputype
mkdir -p $BUILDDIR
cd $BUILDDIR
echo "cmake ../../ $CURDIR -DSHARED=1 -DCMAKE_INSTALL_PREFIX=$1/METIS -DCMAKE_C_COMPILER=$2 -DCMAKE_CXX_COMPILER=$3"
cmake ../../ $CURDIR -DSHARED=1 -DCMAKE_INSTALL_PREFIX=$1/METIS -DCMAKE_C_COMPILER=$2 -DCMAKE_CXX_COMPILER=$3
make -j 4
mkdir $1/METIS
make install


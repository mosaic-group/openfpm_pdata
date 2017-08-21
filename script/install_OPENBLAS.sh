#! /bin/bash

# check if the directory $1/OPENBLAS exist

if [ -d "$1/OPENBLAS" ]; then
  echo "OPENBLAS is already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/OpenBLAS-0.2.19.tar.gz
rm -rf OpenBLAS-0.2.19
tar -xf OpenBLAS-0.2.19.tar.gz
cd OpenBLAS-0.2.19

wget http://ppmcore.mpi-cbg.de/upload/openblas.diff
patch -p1 < openblas.diff

# configuration

make
mkdir $1/OPENBLAS
make install PREFIX=$1/OPENBLAS


# if empty remove the folder
if [ ! "$(ls -A $1/OPENBLAS)" ]; then
   rm -rf $1/OPENBLAS
else
   rm -rf OpenBLAS-0.2.19
   echo 1 > $1/OPENBLAS/version
   exit 0
fi


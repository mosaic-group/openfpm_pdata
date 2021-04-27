#! /bin/bash

source script/discover_os

discover_os

# check if the directory $1/OPENBLAS exist

if [ -d "$1/OPENBLAS" ]; then
  echo "OPENBLAS is already installed"
  exit 0
fi

if [ x"$platform" == x"darwin" ]; then
  rm -rf OpenBLAS-0.3.10
  rm -rf OpenBLAS-0.3.10.tar.gz
  wget http://ppmcore.mpi-cbg.de/upload/OpenBLAS-0.3.10.tar.gz
  tar -xf OpenBLAS-0.3.10.tar.gz
  cd OpenBLAS-0.3.10
else
  rm -rf OpenBLAS-0.3.13
  rm -rf OpenBLAS-0.3.13.tar.gz
  wget http://ppmcore.mpi-cbg.de/upload/OpenBLAS-0.3.13.tar.gz
  tar -xf OpenBLAS-0.3.13.tar.gz
  cd OpenBLAS-0.3.13
fi

#wget http://ppmcore.mpi-cbg.de/upload/openblas.diff
#patch -p1 < openblas.diff

# configuration

make FC=$FC CC=$CC
mkdir $1/OPENBLAS
make install PREFIX=$1/OPENBLAS


# if empty remove the folder
if [ ! "$(ls -A $1/OPENBLAS)" ]; then
   rm -rf $1/OPENBLAS
else
   rm -rf OpenBLAS-0.3.13
   rm -rf OpenBLAS-0.3.10
   echo 3 > $1/OPENBLAS/version
   exit 0
fi

